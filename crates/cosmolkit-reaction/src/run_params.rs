use crate::ReactionCoordinateSelection;
/// Options for ordered multiple-reagent reaction execution.
#[derive(Debug, Clone, PartialEq, Eq, cosmolkit_macros::Configuration)]
pub struct ReactionRunParams {
    /// RDKit default 1000; zero leaves matching and combinations unlimited.
    pub max_products: u32,
    /// Empty selects Auto for each input; otherwise input arity must match.
    #[configuration(none_default)]
    pub coordinate_selections: Vec<ReactionCoordinateSelection>,
    /// CK extension: fill missing product user properties from atom origins.
    /// False preserves native behavior; it never clears existing properties.
    pub copy_atom_properties: bool,
}
impl Default for ReactionRunParams {
    fn default() -> Self {
        Self {
            max_products: 1000,
            coordinate_selections: Vec::new(),
            copy_atom_properties: false,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, cosmolkit_macros::Configuration)]
#[configuration(scalar)]
/// Coordinate selection for executing one reaction template against one input.
pub struct ReactionSingleRunParams {
    #[configuration(none_default)]
    /// Automatic or explicit input-conformer selection.
    pub coordinate_selection: ReactionCoordinateSelection,
}
