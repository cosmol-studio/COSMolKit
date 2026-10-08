use cosmolkit_model::{AtomId, BondId};
#[derive(Debug, thiserror::Error)]
pub enum ReactionProductError {
    #[error(transparent)]
    Coordinate(#[from] cosmolkit_model::CoordinateValidationError),
    #[error(transparent)]
    TopologyEdit(#[from] cosmolkit_model::TopologyEditError),
    #[error(transparent)]
    Adjacency(#[from] cosmolkit_model::AdjacencyError),
    #[error(transparent)]
    RingFinding(#[from] cosmolkit_core::RingFindingError),
    #[error(transparent)]
    SourceUInt(#[from] cosmolkit_core::PropertyUIntReadError),
    #[error(transparent)]
    Valence(#[from] cosmolkit_core::ValenceError),
    #[error(transparent)]
    StereoOrder(#[from] cosmolkit_core::StereoOrderError),
    #[error(transparent)]
    DoubleBondStereo(#[from] cosmolkit_core::DoubleBondStereoError),
    #[error("atom {atom} stereo getter: {source}")]
    StereoGetter {
        atom: AtomId,
        #[source]
        source: cosmolkit_search::SubstructMatchError,
    },
    #[error(transparent)]
    CoordinateSelection(#[from] cosmolkit_smiles::SmilesParseError),
    #[error(transparent)]
    CarrierIdentity(#[from] cosmolkit_model::QueryAtomConversionError),
    #[error("atom {atom}, property {key}: {source}")]
    PropertyInt {
        atom: AtomId,
        key: &'static str,
        #[source]
        source: cosmolkit_core::PropertyIntReadError,
    },
    #[error("atom {atom}, property {key}: {source}")]
    PropertyUInt {
        atom: AtomId,
        key: &'static str,
        #[source]
        source: cosmolkit_core::PropertyUIntReadError,
    },
    #[error("atom {atom} lacks source-required property {key}")]
    MissingProperty { atom: AtomId, key: &'static str },
    #[error(transparent)]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error(transparent)]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error(transparent)]
    BondValue(#[from] cosmolkit_model::BondValueError),
    #[error(transparent)]
    TemplateProperty(#[from] crate::ReactionValidationError),
    #[error(
        "{stage}: {detail} (reactant atom={reactant_atom:?}, product atom={product_atom:?}, bond={bond:?})"
    )]
    Invariant {
        stage: &'static str,
        detail: &'static str,
        reactant_atom: Option<usize>,
        product_atom: Option<usize>,
        bond: Option<BondId>,
    },
    #[error("{kind} row {index} does not fit source unsigned int")]
    RowOverflow { kind: &'static str, index: usize },
    #[error(
        "query-bearing reagent bond {bond} requires the independently excluded query-reactant capability"
    )]
    QueryReactantBond { bond: BondId },
}
