#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum ReactionApplyError {
    #[error(transparent)]
    SourceText(#[from] cosmolkit_core::PropertyStringError),
    #[error(transparent)]
    Initialization(#[from] crate::ReactionInitializationError),
    #[error(
        "only single reactant - single product reactions can be applied; found {reactants} reactants and {products} products"
    )]
    ApplicabilityArity { reactants: usize, products: usize },
    #[error("restricted application cannot add unmapped or new product atom {atom}")]
    AddsProductAtom { atom: cosmolkit_model::AtomId },
    #[error(transparent)]
    Matching(#[from] crate::ReactionRunError),
    #[error(transparent)]
    Product(#[from] crate::ReactionProductError),
    #[error(transparent)]
    Edit(#[from] cosmolkit_model::TopologyEditError),
}
