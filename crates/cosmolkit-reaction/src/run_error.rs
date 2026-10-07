use crate::ReactionInitializationError;

#[derive(Debug, thiserror::Error)]
pub enum ReactionRunError {
    #[error("product set {set}, product template {template}: {source}")]
    Product {
        set: usize,
        template: usize,
        #[source]
        source: crate::ReactionProductError,
    },
    #[error("received {actual} coordinate selectors for {expected} reactants")]
    CoordinateSelectionArity { expected: usize, actual: usize },
    #[error(transparent)]
    Initialization(#[from] ReactionInitializationError),
    #[error("reaction match requires initialized templates")]
    NeedsInitialization,
    #[error("reaction has {expected} reactant templates but received {actual} reactants")]
    ReactantArity { expected: usize, actual: usize },
    #[error("reactant template {index} is outside {count} templates")]
    ReactantTemplateIndex { index: usize, count: usize },
    #[error("reactant {reactant}, template {template} matching: {source}")]
    Matching {
        reactant: usize,
        template: usize,
        #[source]
        source: cosmolkit_search::SubstructMatchError,
    },
    #[error("reactant combination recursion has no levels")]
    EmptyCombinationLevels,
}
