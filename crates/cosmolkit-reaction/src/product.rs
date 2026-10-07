use cosmolkit_core::{RingInfo, ValenceAssignment};
use cosmolkit_model::{AtomId, BondId, CoordinateBlock, MoleculeProperties, TopologyBlock};

/// Authorized borrowed input. This value has no live owner or write authority.
#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct ReactionInput<'a> {
    pub topology: &'a TopologyBlock,
    pub coordinates: &'a CoordinateBlock,
    pub properties: &'a MoleculeProperties,
    /// Current validity-gated facts; absence is checked only by reached getters.
    pub rings: Option<&'a RingInfo>,
    pub valence: Option<&'a ValenceAssignment>,
}

/// One source row of a reconstructed destination, retaining its input index.
/// Duplicate origins represent source-defined one-to-many product mapping.
#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReactionRowOrigin<Row> {
    pub input: usize,
    pub row: Row,
}

/// Complete detached product and its source-derived cache facts.
/// Row origins align with destination row positions; None identifies new rows.
/// The existing runtime validates every source boundary before constructing
/// any live output. This is intentionally distinct from inverse mappings.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionProduct {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub atom_origins: Vec<Option<ReactionRowOrigin<AtomId>>>,
    pub bond_origins: Vec<Option<ReactionRowOrigin<BondId>>>,
    pub valence: ValenceAssignment,
    pub rings: Option<RingInfo>,
}
