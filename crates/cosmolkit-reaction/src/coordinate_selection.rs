/// Coordinate selection for a reaction input or template.
///
/// Auto omits absent coordinates and selects the first source conformer using
/// the canonical source insertion order. Explicit dimension/id selections are
/// the approved detached projection; selection is reached only for source
/// conformers during product construction.
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
