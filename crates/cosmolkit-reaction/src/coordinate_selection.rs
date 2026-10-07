/// Coordinate selection for a reaction input or template.
///
/// Auto omits absent coordinates, selects one stored set, and reports ambiguity
/// for multiple sets. ROOT-approved D4 differs from RDKit's first insertion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ReactionCoordinateSelection {
    #[default]
    Auto,
    TwoD {
        id: usize,
    },
    ThreeD {
        id: usize,
    },
}
impl From<ReactionCoordinateSelection> for cosmolkit_smiles::CxCoordinateSelection {
    fn from(value: ReactionCoordinateSelection) -> Self {
        match value {
            ReactionCoordinateSelection::Auto => Self::Auto,
            ReactionCoordinateSelection::TwoD { id } => Self::TwoD { id },
            ReactionCoordinateSelection::ThreeD { id } => Self::ThreeD { id },
        }
    }
}
