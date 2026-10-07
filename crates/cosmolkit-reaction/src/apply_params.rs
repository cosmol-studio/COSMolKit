#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactionApplyParams {
    pub remove_unmatched_atoms: bool,
}
impl Default for ReactionApplyParams {
    fn default() -> Self {
        Self {
            remove_unmatched_atoms: true,
        }
    }
}
