use crate::{ReactionModelError, ReactionRole};
use cosmolkit_model::{PropertyStore, QueryGraph};
use cosmolkit_search::SubstructMatchParams;

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
    pub(crate) properties: PropertyStore,
}

impl Clone for Reaction {
    fn clone(&self) -> Self {
        // RDKit❗✔️:   ChemicalReaction(const ChemicalReaction &other) : RDProps() { copy(other); }
        let mut result = Self::new();
        result.copy_from_source(self);
        result
    }

    fn clone_from(&mut self, other: &Self) {
        // RDKit❗✔️:   ChemicalReaction &operator=(const ChemicalReaction &other) {
        // RDKit❗✔️:     if (this != &other) {
        // RDKit❗✔️:       copy(other);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     return *this;
        // RDKit❗✔️:   }
        // Safe Rust exclusive/shared references cannot self-alias here.
        self.copy_from_source(other);
    }
}

impl Default for Reaction {
    fn default() -> Self {
        Self::new()
    }
}

impl Reaction {
    fn copy_from_source(&mut self, other: &Self) {
        // RDKit❗✔️:   void copy(const ChemicalReaction &other) {
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
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RDProps &operator=(const RDProps &rhs) {
        // RDKit❗✔️:     if (this == &rhs) {
        // RDKit❗✔️:       return *this;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     d_props = rhs.d_props;
        // RDKit❗✔️:     return *this;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   Dict &operator=(const Dict &other) {
        // RDKit❗✔️:     if (this == &other) {
        // RDKit❗✔️:       return *this;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (_hasNonPodData) {
        // RDKit❗✔️:       reset();
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     if (other._hasNonPodData) {
        // RDKit❗✔️:       std::vector<Pair> data(other._data.size());
        // RDKit❗✔️:       _data.swap(data);
        // RDKit❗✔️:       for (size_t i = 0; i < _data.size(); ++i) {
        // RDKit❗✔️:         _data[i].key = other._data[i].key;
        // RDKit❗✔️:         copy_rdvalue(_data[i].val, other._data[i].val);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     } else {
        // RDKit❗✔️:       _data = other._data;
        // RDKit❗✔️:     }
        // RDKit❗✔️:     _hasNonPodData = other._hasNonPodData;
        // RDKit❗✔️:     return *this;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   RWMol(const ROMol &other, bool quickCopy = false, int confId = -1)
        // RDKit❗✔️:       : ROMol(other, quickCopy, confId) {}
        // RDKit❗✔️:   ROMol(const ROMol &other, bool quickCopy = false, int confId = -1)
        // RDKit❗✔️:       : RDProps() {
        // RDKit❗✔️:     dp_ringInfo = nullptr;
        // RDKit❗✔️:     initFromOther(other, quickCopy, confId);
        // RDKit❗✔️:     numBonds = rdcast<unsigned int>(boost::num_edges(d_graph));
        // RDKit❗✔️:   }
        // Canonical MODEL dictionary/QueryGraph clones own all represented
        // tagged byte properties and default quickCopy=false,confId=-1 deep
        // template copies. No property conversion/filter, query parse, chemical
        // validation, init, template sorting or clone-all-then-clear fallback.
        // Source assignment order is properties, flags, each cleared role in
        // encounter order, then complete parameter assignment. Existing role
        // buffers retain capacity exactly across clear/reserve/repopulate.
        // Cost: O(total represented template/query/property/coordinate data),
        // plus scalar configuration/vector copies. MODEL ordered dictionary
        // indexes and adjacency are the existing transport storage; no second
        // property implementation, mutable native pointer or runtime owner.
        self.properties.clone_from(&other.properties);
        self.needs_init = other.needs_init;
        self.implicit_properties = other.implicit_properties;
        self.reactants.clear();
        self.reactants.reserve(other.reactants.len());
        for template in &other.reactants {
            self.reactants.push(template.clone());
        }
        self.products.clear();
        self.products.reserve(other.products.len());
        for template in &other.products {
            self.products.push(template.clone());
        }
        self.agents.clear();
        self.agents.reserve(other.agents.len());
        for template in &other.agents {
            self.agents.push(template.clone());
        }
        self.match_params.clone_from(&other.match_params);
    }

    pub fn new() -> Self {
        // RDKit❗✔️:   ChemicalReaction() : RDProps() {}
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
            properties: PropertyStore::new(),
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
        // RDKit❗✔️:   const MOL_SPTR_VECT &getReactants() const {
        // RDKit❗✔️:     return this->m_reactantTemplates;
        // RDKit❗✔️:   }
        // Borrow the actual stored sequence: no allocation, validation,
        // initialization, sorting or graph clone; O(1), as native reference.
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
        // RDKit❗✔️:   unsigned int getNumReactantTemplates() const {
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_reactantTemplates.size());
        // RDKit❗✔️:   }
        // Native release rdcast first narrows to unsigned int. The existing
        // Rust projection then widens that count, preserving its low 32 bits.
        (self.reactants.len() as u32) as usize
    }

    pub(crate) fn add_reactant_template_source(&mut self, template: QueryGraph) -> u32 {
        // RDKit❗✔️:   unsigned int addReactantTemplate(ROMOL_SPTR mol) {
        // RDKit❗✔️:     this->df_needsInit = true;
        // RDKit❗✔️:     this->m_reactantTemplates.push_back(mol);
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_reactantTemplates.size());
        // RDKit❗✔️:   }
        // Release rdcast is static_cast: retain the low uint32 count bits.
        // Move the actual detached template without cloning/validating it,
        // after the exact needs-init write. O(1) amortized append, O(R) only
        // when the vector reallocates, matching native vector::push_back.
        self.needs_init = true;
        self.reactants.push(template);
        self.reactants.len() as u32
    }

    pub fn with_reactant_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        // Existing checked immutable value projection owns structural checks
        // and the full reaction copy; the source mutation below is shared.
        validate_template(&template, ReactionRole::Reactant, self.reactants.len())?;
        let mut result = self.clone();
        result.add_reactant_template_source(template);
        Ok(result)
    }

    pub fn product_templates(&self) -> &[QueryGraph] {
        // RDKit❗✔️:   const MOL_SPTR_VECT &getProducts() const { return this->m_productTemplates; }
        // Borrow the actual stored sequence: no allocation, validation,
        // initialization, sorting or graph clone; O(1), as native reference.
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
        // RDKit❗✔️:   unsigned int getNumProductTemplates() const {
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_productTemplates.size());
        // RDKit❗✔️:   }
        // Native release rdcast first narrows to unsigned int. The existing
        // Rust projection then widens that count, preserving its low 32 bits.
        (self.products.len() as u32) as usize
    }

    pub(crate) fn add_product_template_source(&mut self, template: QueryGraph) -> u32 {
        // RDKit❗✔️:   unsigned int addProductTemplate(ROMOL_SPTR mol) {
        // RDKit❗✔️:     this->m_productTemplates.push_back(mol);
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_productTemplates.size());
        // RDKit❗✔️:   }
        // Source does not alter needs-init for this role. Move the actual
        // template, then apply release rdcast's uint32 truncation to size.
        // O(1) amortized, O(n) only on vector growth; no graph deep copy.
        self.products.push(template);
        self.products.len() as u32
    }

    pub fn with_product_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        // The checked immutable projection owns validation and reaction copy.
        validate_template(&template, ReactionRole::Product, self.products.len())?;
        let mut result = self.clone();
        result.add_product_template_source(template);
        Ok(result)
    }

    pub fn agent_templates(&self) -> &[QueryGraph] {
        // RDKit❗✔️:   const MOL_SPTR_VECT &getAgents() const { return this->m_agentTemplates; }
        // Borrow the actual stored sequence: no allocation, validation,
        // initialization, sorting or graph clone; O(1), as native reference.
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
        // RDKit❗✔️:   unsigned int getNumAgentTemplates() const {
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_agentTemplates.size());
        // RDKit❗✔️:   }
        // Native release rdcast first narrows to unsigned int. The existing
        // Rust projection then widens that count, preserving its low 32 bits.
        (self.agents.len() as u32) as usize
    }

    pub(crate) fn add_agent_template_source(&mut self, template: QueryGraph) -> u32 {
        // RDKit❗✔️:   unsigned int addAgentTemplate(ROMOL_SPTR mol) {
        // RDKit❗✔️:     this->m_agentTemplates.push_back(mol);
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_agentTemplates.size());
        // RDKit❗✔️:   }
        // Source does not alter needs-init for this role. Move the actual
        // template, then apply release rdcast's uint32 truncation to size.
        // O(1) amortized, O(n) only on vector growth; no graph deep copy.
        self.agents.push(template);
        self.agents.len() as u32
    }

    pub fn with_agent_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        // The checked immutable projection owns validation and reaction copy.
        validate_template(&template, ReactionRole::Agent, self.agents.len())?;
        let mut result = self.clone();
        result.add_agent_template_source(template);
        Ok(result)
    }

    fn begin_reactant_templates_source(&self) -> (&[QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator beginReactantTemplates() const {
        // RDKit❗✔️:     return this->m_reactantTemplates.begin();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = 0;
        (&self.reactants, offset)
    }

    fn end_reactant_templates_source(&self) -> (&[QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator endReactantTemplates() const {
        // RDKit❗✔️:     return this->m_reactantTemplates.end();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = self.reactants.len();
        (&self.reactants, offset)
    }

    fn begin_product_templates_source(&self) -> (&[QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator beginProductTemplates() const {
        // RDKit❗✔️:     return this->m_productTemplates.begin();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = 0;
        (&self.products, offset)
    }

    fn end_product_templates_source(&self) -> (&[QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator endProductTemplates() const {
        // RDKit❗✔️:     return this->m_productTemplates.end();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = self.products.len();
        (&self.products, offset)
    }

    fn begin_agent_templates_source(&self) -> (&[QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator beginAgentTemplates() const {
        // RDKit❗✔️:     return this->m_agentTemplates.begin();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = 0;
        (&self.agents, offset)
    }

    fn end_agent_templates_source(&self) -> (&[QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator endAgentTemplates() const {
        // RDKit❗✔️:     return this->m_agentTemplates.end();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = self.agents.len();
        (&self.agents, offset)
    }

    fn begin_reactant_templates_mut_source(&mut self) -> (&mut [QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::iterator beginReactantTemplates() {
        // RDKit❗✔️:     return this->m_reactantTemplates.begin();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = 0;
        (&mut self.reactants, offset)
    }

    fn end_reactant_templates_mut_source(&mut self) -> (&mut [QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::iterator endReactantTemplates() {
        // RDKit❗✔️:     return this->m_reactantTemplates.end();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = self.reactants.len();
        (&mut self.reactants, offset)
    }

    fn begin_product_templates_mut_source(&mut self) -> (&mut [QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::iterator beginProductTemplates() {
        // RDKit❗✔️:     return this->m_productTemplates.begin();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = 0;
        (&mut self.products, offset)
    }

    fn end_product_templates_mut_source(&mut self) -> (&mut [QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::iterator endProductTemplates() {
        // RDKit❗✔️:     return this->m_productTemplates.end();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = self.products.len();
        (&mut self.products, offset)
    }

    fn begin_agent_templates_mut_source(&mut self) -> (&mut [QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::iterator beginAgentTemplates() {
        // RDKit❗✔️:     return this->m_agentTemplates.begin();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = 0;
        (&mut self.agents, offset)
    }

    fn end_agent_templates_mut_source(&mut self) -> (&mut [QueryGraph], usize) {
        // RDKit❗✔️:   MOL_SPTR_VECT::iterator endAgentTemplates() {
        // RDKit❗✔️:     return this->m_agentTemplates.end();
        // RDKit❗✔️:   }
        // An iterator position is represented by the actual borrowed sequence
        // plus its source offset. Retaining the whole sequence also preserves
        // decrement/offset access from end, unlike an empty-tail iterator.
        // No template copy, validation, initialization or allocation: O(1).
        let offset = self.agents.len();
        (&mut self.agents, offset)
    }

    pub fn is_initialized(&self) -> bool {
        // RDKit❗✔️:   bool isInitialized() const { return !df_needsInit; }
        !self.needs_init
    }

    pub fn implicit_properties(&self) -> bool {
        // RDKit❗✔️:   bool getImplicitPropertiesFlag() const { return df_implicitProperties; }
        self.implicit_properties
    }

    pub(crate) fn set_implicit_properties_source(&mut self, enabled: bool) {
        // RDKit❗✔️:   void setImplicitPropertiesFlag(bool val) { df_implicitProperties = val; }
        self.implicit_properties = enabled;
    }

    pub fn with_implicit_properties(&self, enabled: bool) -> Self {
        // The immutable value projection separately owns the full copy.
        let mut result = self.clone();
        result.set_implicit_properties_source(enabled);
        result
    }

    pub fn match_params(&self) -> &SubstructMatchParams {
        // RDKit❗✔️:   const SubstructMatchParameters &getSubstructParams() const {
        // RDKit❗✔️:     return d_substructParams;
        // RDKit❗✔️:   }
        &self.match_params
    }

    fn match_params_mut_source(&mut self) -> &mut SubstructMatchParams {
        // RDKit❗✔️:   SubstructMatchParameters &getSubstructParams() { return d_substructParams; }
        &mut self.match_params
    }

    pub fn with_match_params(
        &self,
        params: &SubstructMatchParams,
    ) -> Result<Self, ReactionModelError> {
        // Existing immutable value projection owns copy/assignment, while the
        // native-shaped mutable getter owns direct configuration access.
        let mut result = self.clone();
        result.match_params_mut_source().clone_from(params);
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

#[cfg(test)]
mod complete_reaction_copy_source_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, Conformer2D, PropertyText, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;
    fn template(label: &str) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [
                (PropertyText::from("label"), PropertyValue::from(label)),
                (PropertyText::from("calc"), PropertyValue::Int(7)),
                (
                    PropertyText::from("__computedProps"),
                    PropertyValue::StringVector(vec![PropertyText::from("calc")]),
                ),
            ],
            vec![Conformer2D::new(77, vec![[1.0, 2.0]])],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn reaction() -> Reaction {
        Reaction::from_templates(
            vec![template("r0"), template("r1")],
            vec![template("p0")],
            vec![template("a0"), template("a1")],
        )
        .unwrap()
    }
    #[test]
    fn copies_flags_all_roles_and_parameters_without_initializing_or_validating_chemistry() {
        let mut original = reaction();
        original.needs_init = false;
        original.implicit_properties = true;
        original.match_params.max_matches = 17;
        original.match_params.num_threads = -2;
        original.match_params.use_chirality = true;
        original.match_params.atom_properties = vec!["z".into(), "a".into()];
        let copied = original.clone();
        assert!(!copied.needs_init);
        assert!(copied.implicit_properties);
        assert_eq!(copied.reactants, original.reactants);
        assert_eq!(copied.products, original.products);
        assert_eq!(copied.agents, original.agents);
        assert_eq!(
            format!("{:?}", copied.match_params),
            format!("{:?}", original.match_params)
        );
    }
    #[test]
    fn typed_ordered_properties_preserve_empty_nul_non_utf8_keys_and_computed_metadata() {
        let mut original = reaction();
        let raw = PropertyText::from(vec![b'x', 0, 0xff]);
        original.properties = PropertyStore::from_records([
            (PropertyText::from("z"), PropertyValue::UInt(u32::MAX)),
            (PropertyText::from(""), PropertyValue::Bool(false)),
            (
                raw.clone(),
                PropertyValue::String(PropertyText::from(vec![0xff, 0, b'v'])),
            ),
            (
                PropertyText::from("numbers"),
                PropertyValue::IntVector(vec![i32::MIN, 0, i32::MAX]),
            ),
            (
                PropertyText::from("__computedProps"),
                PropertyValue::StringVector(vec![raw.clone(), PropertyText::from("z")]),
            ),
        ]);
        let mut copied = original.clone();
        assert_eq!(copied.properties, original.properties);
        assert_eq!(
            copied
                .properties
                .ordered()
                .map(|(key, _)| key.clone())
                .collect::<Vec<_>>(),
            vec![
                "z".into(),
                "".into(),
                raw.clone(),
                "numbers".into(),
                "__computedProps".into()
            ]
        );
        assert_eq!(
            copied.properties.get(raw.as_bytes()),
            original.properties.get(raw.as_bytes())
        );
        copied.properties = PropertyStore::new();
        assert!(copied.properties.ordered().next().is_none());
        assert_eq!(original.properties.ordered().len(), 5);
    }
    #[test]
    fn assignment_replaces_existing_roles_and_properties_but_reuses_each_role_buffer() {
        let source = reaction();
        let mut target = reaction();
        target.reactants.reserve(40);
        target.products.reserve(40);
        target.agents.reserve(40);
        let pointers = (
            target.reactants.as_ptr(),
            target.products.as_ptr(),
            target.agents.as_ptr(),
        );
        target.properties = PropertyStore::from_records([("stale".into(), PropertyValue::Int(1))]);
        target.clone_from(&source);
        assert_eq!(
            (
                target.reactants.as_ptr(),
                target.products.as_ptr(),
                target.agents.as_ptr()
            ),
            pointers
        );
        assert_eq!(target.reactants, source.reactants);
        assert_eq!(target.products, source.products);
        assert_eq!(target.agents, source.agents);
        assert!(target.properties.get(b"stale").is_none());
        assert_eq!(target.properties, source.properties);
    }
    #[test]
    fn deep_template_copies_retain_coordinates_and_computed_props_without_shared_mutable_values() {
        let original = reaction();
        let mut copied = original.clone();
        assert_ne!(copied.reactants.as_ptr(), original.reactants.as_ptr());
        assert_ne!(
            copied.reactants[0].atoms().as_ptr(),
            original.reactants[0].atoms().as_ptr()
        );
        assert_eq!(copied.reactants[0].conformers_2d()[0].id(), 77);
        assert_eq!(
            copied.reactants[0].conformers_2d()[0].coordinates(),
            &[[1., 2.]]
        );
        assert_eq!(
            copied.reactants[0].prop("calc"),
            Some(&PropertyValue::Int(7))
        );
        copied.reactants[0]
            .set_prop("calc", PropertyValue::Int(99))
            .unwrap();
        assert_eq!(
            original.reactants[0].prop("calc"),
            Some(&PropertyValue::Int(7))
        );
        assert_eq!(
            original.reactants[0].prop("__computedProps"),
            Some(&PropertyValue::StringVector(vec!["calc".into()]))
        );
    }
    #[test]
    fn empty_source_clears_all_existing_templates_and_overwrites_every_flag_and_parameter() {
        let mut target = reaction();
        target.needs_init = false;
        target.implicit_properties = true;
        target.match_params.max_matches = 1;
        target.clone_from(&Reaction::new());
        assert!(
            target.reactants.is_empty() && target.products.is_empty() && target.agents.is_empty()
        );
        assert!(target.needs_init);
        assert!(!target.implicit_properties);
        assert_eq!(
            format!("{:?}", target.match_params),
            format!("{:?}", SubstructMatchParams::default())
        );
    }
}

#[cfg(test)]
mod complete_add_reactant_template_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, PropertyText, PropertyValue, QueryAtom, QueryGraphError,
    };
    use cosmolkit_types::Element;
    fn template(label: &str) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [(PropertyText::from("label"), PropertyValue::from(label))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn source_append_moves_exact_template_and_sets_init_before_returning_monotonic_counts() {
        let mut reaction =
            Reaction::from_templates(vec![], vec![template("p")], vec![template("a")]).unwrap();
        reaction.implicit_properties = true;
        reaction.match_params.max_matches = 19;
        let product_ptr = reaction.products.as_ptr();
        let agent_ptr = reaction.agents.as_ptr();
        for (i, label) in ["r0", "r1", "r2"].into_iter().enumerate() {
            reaction.needs_init = false;
            let graph = template(label);
            let atoms_ptr = graph.atoms().as_ptr();
            assert_eq!(reaction.add_reactant_template_source(graph), i as u32 + 1);
            assert!(reaction.needs_init);
            assert_eq!(reaction.reactants[i].atoms().as_ptr(), atoms_ptr);
            assert_eq!(
                reaction.reactants[i].prop("label"),
                Some(&PropertyValue::from(label))
            );
        }
        assert_eq!(
            (reaction.products.as_ptr(), reaction.agents.as_ptr()),
            (product_ptr, agent_ptr)
        );
        assert!(reaction.implicit_properties);
        assert_eq!(reaction.match_params.max_matches, 19);
    }
    #[test]
    fn immutable_projection_keeps_source_and_reuses_the_same_template_append_owner() {
        let mut source = Reaction::from_templates(
            vec![template("first")],
            vec![template("p")],
            vec![template("a")],
        )
        .unwrap();
        source.needs_init = false;
        let graph = template("last");
        let atoms_ptr = graph.atoms().as_ptr();
        let result = source.with_reactant_template(graph).unwrap();
        assert!(!source.needs_init);
        assert_eq!(source.reactants.len(), 1);
        assert!(result.needs_init);
        assert_eq!(result.reactants.len(), 2);
        assert_eq!(result.reactants[1].atoms().as_ptr(), atoms_ptr);
        assert_eq!(result.products, source.products);
        assert_eq!(result.agents, source.agents);
    }
    #[test]
    fn structural_validation_stays_in_checked_projection_without_becoming_a_source_precondition() {
        let mut source = Reaction::new();
        source.needs_init = false;
        let mut graph = template("invalid");
        graph.atoms_mut()[0] = QueryAtom::new(AtomId::new(99), AtomSpec::new(Element::C));
        assert!(matches!(
            source.with_reactant_template(graph.clone()),
            Err(ReactionModelError::Template {
                role: ReactionRole::Reactant,
                template: 0,
                source: QueryGraphError::AtomIdMismatch { .. }
            })
        ));
        assert!(!source.needs_init);
        assert!(source.reactants.is_empty());
        assert_eq!(source.add_reactant_template_source(graph), 1);
        assert!(source.needs_init);
        assert!(matches!(
            source.reactants[0].validate(),
            Err(QueryGraphError::AtomIdMismatch { .. })
        ));
    }
}

#[cfg(test)]
mod complete_add_agent_template_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, PropertyText, PropertyValue, QueryAtom, QueryGraphError,
    };
    use cosmolkit_types::Element;
    fn template(label: &str) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [(PropertyText::from("label"), PropertyValue::from(label))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn source_append_moves_template_returns_counts_and_preserves_both_init_states() {
        let mut reaction = Reaction::from_templates(
            vec![template("r")],
            vec![template("p")],
            vec![template("a")],
        )
        .unwrap();
        reaction.implicit_properties = true;
        reaction.match_params.max_matches = 19;
        let reactant_ptr = reaction.reactants.as_ptr();
        let other_ptr = reaction.products.as_ptr();
        for (i, needs_init) in [false, true, false].into_iter().enumerate() {
            reaction.needs_init = needs_init;
            let graph = template("new");
            let atoms_ptr = graph.atoms().as_ptr();
            assert_eq!(reaction.add_agent_template_source(graph), i as u32 + 2);
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.agents[i + 1].atoms().as_ptr(), atoms_ptr);
            assert_eq!(
                reaction.agents[i + 1].prop("label"),
                Some(&PropertyValue::from("new"))
            );
        }
        assert_eq!(
            (reaction.reactants.as_ptr(), reaction.products.as_ptr()),
            (reactant_ptr, other_ptr)
        );
        assert!(reaction.implicit_properties);
        assert_eq!(reaction.match_params.max_matches, 19);
    }
    #[test]
    fn immutable_projection_keeps_original_and_initialized_state_without_reordering() {
        let mut source = Reaction::from_templates(
            vec![template("r")],
            vec![template("p")],
            vec![template("a")],
        )
        .unwrap();
        source.needs_init = false;
        let graph = template("last");
        let atoms_ptr = graph.atoms().as_ptr();
        let result = source.with_agent_template(graph).unwrap();
        assert!(!source.needs_init);
        assert_eq!(source.agents.len(), 1);
        assert!(!result.needs_init);
        assert_eq!(result.agents.len(), 2);
        assert_eq!(result.agents[1].atoms().as_ptr(), atoms_ptr);
        assert_eq!(result.reactants, source.reactants);
        assert_eq!(result.products, source.products);
        assert_eq!(result.agents[0], source.agents[0]);
    }
    #[test]
    fn checked_projection_validates_without_adding_source_preconditions() {
        let mut source = Reaction::new();
        source.needs_init = false;
        let mut graph = template("invalid");
        graph.atoms_mut()[0] = QueryAtom::new(AtomId::new(99), AtomSpec::new(Element::C));
        assert!(matches!(
            source.with_agent_template(graph.clone()),
            Err(ReactionModelError::Template {
                role: ReactionRole::Agent,
                template: 0,
                source: QueryGraphError::AtomIdMismatch { .. }
            })
        ));
        assert!(!source.needs_init);
        assert!(source.agents.is_empty());
        assert_eq!(source.add_agent_template_source(graph), 1);
        assert!(!source.needs_init);
        assert!(matches!(
            source.agents[0].validate(),
            Err(QueryGraphError::AtomIdMismatch { .. })
        ));
    }
}

#[cfg(test)]
mod complete_add_product_template_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, PropertyText, PropertyValue, QueryAtom, QueryGraphError,
    };
    use cosmolkit_types::Element;
    fn template(label: &str) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [(PropertyText::from("label"), PropertyValue::from(label))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn source_append_moves_template_returns_counts_and_preserves_both_init_states() {
        let mut reaction = Reaction::from_templates(
            vec![template("r")],
            vec![template("p")],
            vec![template("a")],
        )
        .unwrap();
        reaction.implicit_properties = true;
        reaction.match_params.max_matches = 19;
        let reactant_ptr = reaction.reactants.as_ptr();
        let other_ptr = reaction.agents.as_ptr();
        for (i, needs_init) in [false, true, false].into_iter().enumerate() {
            reaction.needs_init = needs_init;
            let graph = template("new");
            let atoms_ptr = graph.atoms().as_ptr();
            assert_eq!(reaction.add_product_template_source(graph), i as u32 + 2);
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.products[i + 1].atoms().as_ptr(), atoms_ptr);
            assert_eq!(
                reaction.products[i + 1].prop("label"),
                Some(&PropertyValue::from("new"))
            );
        }
        assert_eq!(
            (reaction.reactants.as_ptr(), reaction.agents.as_ptr()),
            (reactant_ptr, other_ptr)
        );
        assert!(reaction.implicit_properties);
        assert_eq!(reaction.match_params.max_matches, 19);
    }
    #[test]
    fn immutable_projection_keeps_original_and_initialized_state_without_reordering() {
        let mut source = Reaction::from_templates(
            vec![template("r")],
            vec![template("p")],
            vec![template("a")],
        )
        .unwrap();
        source.needs_init = false;
        let graph = template("last");
        let atoms_ptr = graph.atoms().as_ptr();
        let result = source.with_product_template(graph).unwrap();
        assert!(!source.needs_init);
        assert_eq!(source.products.len(), 1);
        assert!(!result.needs_init);
        assert_eq!(result.products.len(), 2);
        assert_eq!(result.products[1].atoms().as_ptr(), atoms_ptr);
        assert_eq!(result.reactants, source.reactants);
        assert_eq!(result.agents, source.agents);
        assert_eq!(result.products[0], source.products[0]);
    }
    #[test]
    fn checked_projection_validates_without_adding_source_preconditions() {
        let mut source = Reaction::new();
        source.needs_init = false;
        let mut graph = template("invalid");
        graph.atoms_mut()[0] = QueryAtom::new(AtomId::new(99), AtomSpec::new(Element::C));
        assert!(matches!(
            source.with_product_template(graph.clone()),
            Err(ReactionModelError::Template {
                role: ReactionRole::Product,
                template: 0,
                source: QueryGraphError::AtomIdMismatch { .. }
            })
        ));
        assert!(!source.needs_init);
        assert!(source.products.is_empty());
        assert_eq!(source.add_product_template_source(graph), 1);
        assert!(!source.needs_init);
        assert!(matches!(
            source.products[0].validate(),
            Err(QueryGraphError::AtomIdMismatch { .. })
        ));
    }
}

#[cfg(test)]
mod complete_template_sequence_access_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn sequence_access_borrows_actual_storage_in_order_without_validating_or_initializing() {
        let mut reaction = Reaction::new();
        assert!(reaction.reactant_templates().is_empty());
        assert!(reaction.agent_templates().is_empty());
        assert!(reaction.product_templates().is_empty());
        for (role, count) in [
            (ReactionRole::Reactant, 3),
            (ReactionRole::Product, 2),
            (ReactionRole::Agent, 1),
        ] {
            let rows = match role {
                ReactionRole::Reactant => &mut reaction.reactants,
                ReactionRole::Product => &mut reaction.products,
                ReactionRole::Agent => &mut reaction.agents,
            };
            for i in 0..count {
                let mut graph = QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap();
                graph.atoms_mut()[0] =
                    QueryAtom::new(AtomId::new(90 + i), AtomSpec::new(Element::C));
                rows.push(graph);
            }
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            for (actual, stored) in [
                (reaction.reactant_templates(), &reaction.reactants),
                (reaction.product_templates(), &reaction.products),
                (reaction.agent_templates(), &reaction.agents),
            ] {
                assert_eq!(actual.as_ptr(), stored.as_ptr());
                assert_eq!(actual.len(), stored.len());
                for (i, graph) in actual.iter().enumerate() {
                    assert_eq!(graph.atoms()[0].id().index(), 90 + i);
                    assert!(graph.validate().is_err());
                }
            }
            assert_eq!(reaction.needs_init, needs_init);
        }
    }
}

#[cfg(test)]
mod complete_begin_reactant_templates_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.begin_reactant_templates_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.reactants.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.reactants.as_ptr();
            let atom_ptr = reaction.reactants[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.begin_reactant_templates_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 0);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                assert_eq!(rows[0].atoms()[0].id().index(), 0);
            }
            assert_eq!(reaction.needs_init, needs_init);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_end_reactant_templates_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.end_reactant_templates_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.reactants.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.reactants.as_ptr();
            let atom_ptr = reaction.reactants[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.end_reactant_templates_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 3);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                assert_eq!(rows[0].atoms()[0].id().index(), 0);
            }
            assert_eq!(reaction.needs_init, needs_init);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_begin_product_templates_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.begin_product_templates_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.products.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.products.as_ptr();
            let atom_ptr = reaction.products[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.begin_product_templates_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 0);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                assert_eq!(rows[0].atoms()[0].id().index(), 0);
            }
            assert_eq!(reaction.needs_init, needs_init);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_end_product_templates_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.end_product_templates_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.products.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.products.as_ptr();
            let atom_ptr = reaction.products[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.end_product_templates_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 3);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                assert_eq!(rows[0].atoms()[0].id().index(), 0);
            }
            assert_eq!(reaction.needs_init, needs_init);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_begin_agent_templates_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.begin_agent_templates_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.agents.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.agents.as_ptr();
            let atom_ptr = reaction.agents[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.begin_agent_templates_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 0);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                assert_eq!(rows[0].atoms()[0].id().index(), 0);
            }
            assert_eq!(reaction.needs_init, needs_init);
        }
        assert!(reaction.products.is_empty());
    }
}

#[cfg(test)]
mod complete_end_agent_templates_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.end_agent_templates_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.agents.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.agents.as_ptr();
            let atom_ptr = reaction.agents[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.end_agent_templates_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 3);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                assert_eq!(rows[0].atoms()[0].id().index(), 0);
            }
            assert_eq!(reaction.needs_init, needs_init);
        }
        assert!(reaction.products.is_empty());
    }
}

#[cfg(test)]
mod complete_begin_reactant_templates_mut_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.begin_reactant_templates_mut_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.reactants.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.reactants.as_ptr();
            let atom_ptr = reaction.reactants[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.begin_reactant_templates_mut_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 0);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                rows[0].atoms_mut()[0] = QueryAtom::new(AtomId::new(42), AtomSpec::new(Element::N));
            }
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.reactants[0].atoms()[0].id().index(), 42);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_end_reactant_templates_mut_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.end_reactant_templates_mut_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.reactants.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.reactants.as_ptr();
            let atom_ptr = reaction.reactants[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.end_reactant_templates_mut_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 3);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                rows[0].atoms_mut()[0] = QueryAtom::new(AtomId::new(42), AtomSpec::new(Element::N));
            }
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.reactants[0].atoms()[0].id().index(), 42);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_begin_product_templates_mut_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.begin_product_templates_mut_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.products.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.products.as_ptr();
            let atom_ptr = reaction.products[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.begin_product_templates_mut_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 0);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                rows[0].atoms_mut()[0] = QueryAtom::new(AtomId::new(42), AtomSpec::new(Element::N));
            }
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.products[0].atoms()[0].id().index(), 42);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_end_product_templates_mut_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.end_product_templates_mut_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.products.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.products.as_ptr();
            let atom_ptr = reaction.products[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.end_product_templates_mut_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 3);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                rows[0].atoms_mut()[0] = QueryAtom::new(AtomId::new(42), AtomSpec::new(Element::N));
            }
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.products[0].atoms()[0].id().index(), 42);
        }
        assert!(reaction.agents.is_empty());
    }
}

#[cfg(test)]
mod complete_begin_agent_templates_mut_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.begin_agent_templates_mut_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.agents.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.agents.as_ptr();
            let atom_ptr = reaction.agents[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.begin_agent_templates_mut_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 0);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                rows[0].atoms_mut()[0] = QueryAtom::new(AtomId::new(42), AtomSpec::new(Element::N));
            }
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.agents[0].atoms()[0].id().index(), 42);
        }
        assert!(reaction.products.is_empty());
    }
}

#[cfg(test)]
mod complete_end_agent_templates_mut_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    #[test]
    fn cursor_preserves_storage_position_order_and_init_without_copy_or_validation() {
        let mut reaction = Reaction::new();
        {
            let (rows, offset) = reaction.end_agent_templates_mut_source();
            assert_eq!(offset, 0);
            assert!(rows.is_empty());
        }
        for element in [Element::C, Element::O, Element::N] {
            reaction.agents.push(
                QueryGraph::from_parts(
                    vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(element))],
                    vec![],
                    [],
                    vec![],
                    vec![],
                    vec![],
                )
                .unwrap(),
            );
        }
        for needs_init in [false, true] {
            reaction.needs_init = needs_init;
            let storage_ptr = reaction.agents.as_ptr();
            let atom_ptr = reaction.agents[2].atoms().as_ptr();
            {
                let (rows, offset) = reaction.end_agent_templates_mut_source();
                assert_eq!(rows.as_ptr(), storage_ptr);
                assert_eq!(rows.len(), 3);
                assert_eq!(offset, 3);
                assert_eq!(rows[2].atoms().as_ptr(), atom_ptr);
                // The complete sequence remains accessible from both endpoints.
                assert_eq!(rows[1].atoms()[0].id().index(), 0);
                rows[0].atoms_mut()[0] = QueryAtom::new(AtomId::new(42), AtomSpec::new(Element::N));
            }
            assert_eq!(reaction.needs_init, needs_init);
            assert_eq!(reaction.agents[0].atoms()[0].id().index(), 42);
        }
        assert!(reaction.products.is_empty());
    }
}

#[cfg(test)]
mod complete_template_public_value_queries {
    use super::*;
    fn empty_template() -> QueryGraph {
        QueryGraph::from_parts(vec![], vec![], [], vec![], vec![], vec![]).unwrap()
    }
    #[test]
    fn queries_follow_role_append_copy_and_configuration_without_changing_original() {
        let mut source = Reaction::new();
        assert_eq!(
            (
                source.num_reactant_templates(),
                source.num_product_templates(),
                source.num_agent_templates()
            ),
            (0, 0, 0)
        );
        assert!(!source.is_initialized());
        assert!(!source.implicit_properties());
        source = source.with_reactant_template(empty_template()).unwrap();
        source = source.with_product_template(empty_template()).unwrap();
        source = source.with_agent_template(empty_template()).unwrap();
        source.needs_init = false;
        let configured = source.with_implicit_properties(true);
        assert!(configured.is_initialized());
        assert!(configured.implicit_properties());
        assert!(source.is_initialized());
        assert!(!source.implicit_properties());
        let products = configured.with_product_template(empty_template()).unwrap();
        let agents = products.with_agent_template(empty_template()).unwrap();
        assert!(agents.is_initialized());
        assert!(agents.implicit_properties());
        let reactants = agents.with_reactant_template(empty_template()).unwrap();
        assert!(!reactants.is_initialized());
        assert!(reactants.implicit_properties());
        assert_eq!(
            (
                reactants.num_reactant_templates(),
                reactants.num_product_templates(),
                reactants.num_agent_templates()
            ),
            (2, 2, 2)
        );
        assert_eq!(
            (
                source.num_reactant_templates(),
                source.num_product_templates(),
                source.num_agent_templates()
            ),
            (1, 1, 1)
        );
    }
}

#[cfg(test)]
mod complete_substruct_params_access {
    use super::*;
    #[test]
    fn params_borrow_stored_value_and_immutable_replacement_preserves_original_without_validation()
    {
        let mut source = Reaction::new();
        source.needs_init = false;
        assert!(std::ptr::eq(source.match_params(), &source.match_params));
        let mut params = source.match_params().clone();
        params.max_matches = 0;
        params.num_threads = -3;
        params.recursion_possible = false;
        let replaced = source.with_match_params(&params).unwrap();
        assert_eq!(replaced.match_params().max_matches, 0);
        assert_eq!(replaced.match_params().num_threads, -3);
        assert!(!replaced.match_params().recursion_possible);
        assert_ne!(source.match_params().max_matches, 0);
        assert_eq!(source.match_params().num_threads, 1);
        assert!(source.is_initialized());
        assert!(replaced.is_initialized());
        assert!(std::ptr::eq(
            replaced.match_params(),
            &replaced.match_params
        ));
    }
}
