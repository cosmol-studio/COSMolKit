/// Source removal threshold. The name is retained from RDKit even though
/// the predicate compares mapped-property count divided by heavy-atom count.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactionTemplateRemovalParams {
    pub threshold_unmapped_atoms: f64,
    pub move_to_agent_templates: bool,
}
impl Default for ReactionTemplateRemovalParams {
    fn default() -> Self {
        // RDKit❗✔️: void removeUnmappedReactantTemplates(double thresholdUnmappedAtoms = 0.2,
        // RDKit❗✔️:                                        bool moveToAgentTemplates = true,
        // RDKit❗✔️:                                        MOL_SPTR_VECT *targetVector = nullptr);
        Self {
            threshold_unmapped_atoms: 0.2,
            move_to_agent_templates: true,
        }
    }
}
