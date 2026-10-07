use crate::ReactionCoordinateSelection;
/// Options for ordered multiple-reagent reaction execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReactionRunParams {
    /// RDKit default 1000; zero leaves matching and combinations unlimited.
    pub max_products: u32,
    /// Empty selects Auto for each input; otherwise input arity must match.
    pub coordinate_selections: Vec<ReactionCoordinateSelection>,
}
impl Default for ReactionRunParams {
    fn default() -> Self {
        Self {
            max_products: 1000,
            coordinate_selections: Vec::new(),
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ReactionSingleRunParams {
    pub coordinate_selection: ReactionCoordinateSelection,
}
