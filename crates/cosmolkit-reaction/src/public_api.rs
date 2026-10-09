//! Public value APIs delegate to the existing detached reaction owners.
use crate::*;
use cosmolkit_model::PropertyText;
impl ReactionCoordinateSelection {
    pub fn auto() -> Self {
        Self::Auto
    }
    pub fn two_d(id: usize) -> Self {
        Self::TwoD { id }
    }
    pub fn three_d(id: usize) -> Self {
        Self::ThreeD { id }
    }
    pub fn id(&self) -> Option<usize> {
        match *self {
            Self::Auto => None,
            Self::TwoD { id } | Self::ThreeD { id } => Some(id),
        }
    }
    pub fn is_auto(&self) -> bool {
        matches!(self, Self::Auto)
    }
    pub fn is_2d(&self) -> bool {
        matches!(self, Self::TwoD { .. })
    }
    pub fn is_3d(&self) -> bool {
        matches!(self, Self::ThreeD { .. })
    }
}
impl Reaction {
    pub fn from_smirks(text: &str) -> Result<Self, ReactionParseError> {
        parse_smirks(text)
    }
    pub fn from_smirks_with_params(
        text: &str,
        params: &ReactionParseParams,
    ) -> Result<Self, ReactionParseError> {
        parse_smirks_with_params(text, params)
    }
    pub fn to_smirks(&self) -> Result<PropertyText, ReactionWriteError> {
        self.to_smirks_with_params(&ReactionWriteParams::default())
    }
    pub fn to_smirks_with_params(
        &self,
        params: &ReactionWriteParams,
    ) -> Result<PropertyText, ReactionWriteError> {
        crate::write::write_smirks_with_params(self, params)
    }
    pub fn to_cx_smirks(&self) -> Result<PropertyText, ReactionWriteError> {
        self.to_cx_smirks_with_params(&ReactionWriteParams::default())
    }
    pub fn to_cx_smirks_with_params(
        &self,
        params: &ReactionWriteParams,
    ) -> Result<PropertyText, ReactionWriteError> {
        crate::write::write_cx_smirks_with_params(self, params, params.cx_fields)
    }
    pub fn validate(&self) -> Result<ReactionValidationReport, ReactionValidationError> {
        self.validate_with_params(&ReactionValidationParams::default())
    }
    pub fn validate_with_params(
        &self,
        params: &ReactionValidationParams,
    ) -> Result<ReactionValidationReport, ReactionValidationError> {
        validate_reaction(self, params)
    }
    pub fn with_initialized(&self) -> Result<Self, ReactionInitializationError> {
        self.with_initialized_with_params(&ReactionValidationParams::default())
    }
    pub fn with_initialized_with_params(
        &self,
        params: &ReactionValidationParams,
    ) -> Result<Self, ReactionInitializationError> {
        initialize_reaction(self, params)
    }
    pub fn without_agents(&self) -> ReactionTemplateRemoval {
        without_agents(self)
    }
    pub fn without_unmapped_reactants(&self) -> ReactionTemplateRemoval {
        self.without_unmapped_reactants_with_params(&ReactionTemplateRemovalParams::default())
    }
    pub fn without_unmapped_reactants_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        without_unmapped_reactants(self, params)
    }
    pub fn without_unmapped_products(&self) -> ReactionTemplateRemoval {
        self.without_unmapped_products_with_params(&ReactionTemplateRemovalParams::default())
    }
    pub fn without_unmapped_products_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        without_unmapped_products(self, params)
    }
}
impl ReactionTemplateRemoval {
    pub fn reaction(&self) -> &Reaction {
        &self.reaction
    }
    pub fn removed_templates(&self) -> &[cosmolkit_model::QueryGraph] {
        &self.removed_templates
    }
}
impl ReactionValidationReport {
    pub fn warnings(&self) -> &[ReactionValidationIssue] {
        &self.warnings
    }
    pub fn errors(&self) -> &[ReactionValidationIssue] {
        &self.errors
    }
}
impl ReactionValidationIssue {
    pub fn kind(&self) -> ReactionValidationIssueKind {
        self.kind
    }
    pub fn severity(&self) -> ReactionValidationSeverity {
        self.severity
    }
    pub fn role(&self) -> Option<ReactionRole> {
        self.role
    }
    pub fn template(&self) -> Option<usize> {
        self.template
    }
    pub fn atom(&self) -> Option<cosmolkit_model::AtomId> {
        self.atom
    }
    pub fn map(&self) -> Option<i32> {
        self.map
    }
    pub fn maps(&self) -> &[i32] {
        &self.maps
    }
    pub fn detail(&self) -> &str {
        &self.detail
    }
}
impl ReactionParseParams {
    pub fn new(
        use_smiles: bool,
        sanitize: bool,
        replacements: std::collections::BTreeMap<String, String>,
        allow_cxsmiles: bool,
        strict_cxsmiles: bool,
    ) -> Self {
        Self {
            use_smiles,
            sanitize,
            replacements,
            allow_cxsmiles,
            strict_cxsmiles,
        }
    }
    pub fn use_smiles(&self) -> bool {
        self.use_smiles
    }
    pub fn sanitize(&self) -> bool {
        self.sanitize
    }
    pub fn replacements(&self) -> &std::collections::BTreeMap<String, String> {
        &self.replacements
    }
    pub fn allow_cxsmiles(&self) -> bool {
        self.allow_cxsmiles
    }
    pub fn strict_cxsmiles(&self) -> bool {
        self.strict_cxsmiles
    }
}
impl ReactionValidationParams {
    pub fn new(silent: bool) -> Self {
        Self { silent }
    }
    pub fn silent(&self) -> bool {
        self.silent
    }
}
impl ReactionSingleRunParams {
    pub fn new(coordinate_selection: crate::ReactionCoordinateSelection) -> Self {
        Self {
            coordinate_selection,
        }
    }
    pub fn coordinate_selection(&self) -> crate::ReactionCoordinateSelection {
        self.coordinate_selection
    }
}
impl ReactionRunParams {
    pub fn new(
        max_products: u32,
        coordinate_selections: Vec<crate::ReactionCoordinateSelection>,
        copy_atom_properties: bool,
    ) -> Self {
        Self {
            max_products,
            coordinate_selections,
            copy_atom_properties,
        }
    }
    pub fn max_products(&self) -> u32 {
        self.max_products
    }
    pub fn coordinate_selections(&self) -> &[crate::ReactionCoordinateSelection] {
        &self.coordinate_selections
    }
    pub fn copy_atom_properties(&self) -> bool {
        self.copy_atom_properties
    }
}
impl ReactionApplyParams {
    pub fn new(remove_unmatched_atoms: bool) -> Self {
        Self {
            remove_unmatched_atoms,
        }
    }
    pub fn remove_unmatched_atoms(&self) -> bool {
        self.remove_unmatched_atoms
    }
}
impl ReactionTemplateRemovalParams {
    pub fn new(threshold_unmapped_atoms: f64, move_to_agent_templates: bool) -> Self {
        Self {
            threshold_unmapped_atoms,
            move_to_agent_templates,
        }
    }
    pub fn threshold_unmapped_atoms(&self) -> f64 {
        self.threshold_unmapped_atoms
    }
    pub fn move_to_agent_templates(&self) -> bool {
        self.move_to_agent_templates
    }
}
impl ReactionWriteParams {
    pub fn new(
        canonical: bool,
        isomeric_smiles: bool,
        rooted_at_atom: Option<usize>,
        include_dative_bonds: bool,
        include_cx: bool,
        cx_fields: crate::CxSmilesFields,
        coordinate_selections: Vec<crate::ReactionCoordinateSelection>,
    ) -> Self {
        Self {
            canonical,
            isomeric_smiles,
            rooted_at_atom,
            include_dative_bonds,
            include_cx,
            cx_fields,
            coordinate_selections,
        }
    }
    pub fn canonical(&self) -> bool {
        self.canonical
    }
    pub fn isomeric_smiles(&self) -> bool {
        self.isomeric_smiles
    }
    pub fn rooted_at_atom(&self) -> Option<usize> {
        self.rooted_at_atom
    }
    pub fn include_dative_bonds(&self) -> bool {
        self.include_dative_bonds
    }
    pub fn include_cx(&self) -> bool {
        self.include_cx
    }
    pub fn cx_fields(&self) -> crate::CxSmilesFields {
        self.cx_fields
    }
    pub fn coordinate_selections(&self) -> &[crate::ReactionCoordinateSelection] {
        &self.coordinate_selections
    }
}
