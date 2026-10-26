#[derive(Debug, Clone, Copy, PartialEq, Eq, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
/// Options for applying a reaction to one molecular value.
pub struct ReactionApplyParams {
    /// Remove atoms not matched by the reaction's reactant template.
    pub remove_unmatched_atoms: bool,
}
impl Default for ReactionApplyParams {
    fn default() -> Self {
        Self {
            remove_unmatched_atoms: true,
        }
    }
}
