//! Shared concrete molecular value types for COSMolKit algorithms.
//!
//! This crate owns the data model below a live `Molecule`. Values may be
//! edited while detached; installation into a live molecule remains the
//! responsibility of the parent runtime crate. Generic query payloads on
//! `Atom` and `Bond` are opaque here, so this crate does not depend on search.

mod adjacency;
mod atom;
mod bond;
mod cip;
mod coordinates;
mod properties;
pub mod query;
mod sgroup;
mod stereo_group;
mod topology;

pub use adjacency::{AdjacencyError, AdjacencyList, NeighborRef};
pub use atom::{Atom, AtomPdbResidueInfo, AtomSpec};
pub use bond::{Bond, BondSpec};
pub use cip::{CipDescriptor, CipDescriptorError};
pub use coordinates::{
    Conformer2D, Conformer3D, ConformerStore, CoordinateBlock, CoordinateDimension,
    CoordinateValidationError,
};
pub use properties::{MoleculeProperties, PropertyStore, SdfPropertyList, SdfPropertyListTarget};
pub use query::{
    AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, BondQueryPredicate,
    QueryAtom, QueryBond, QueryGraph, QueryGraphError, QueryNode, RecursiveStructureQuery,
};
pub use sgroup::{
    SGroupAttachPoint, SGroupBondRole, SGroupBracket, SGroupBracketStyle, SGroupCState,
    SGroupConnection, SGroupData, SGroupDisplay, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind,
};
pub use stereo_group::{StereoGroup, StereoGroupKind};
pub use topology::{
    AtomMapping, BondMapping, MappingValidationError, TopologyBlock, TopologyMapping,
    TopologyValidationError,
};

use std::fmt;

pub use cosmolkit_types::{
    BondDirection, BondOrder, BondStereo, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element,
    ElementInfo, ElementParseError, Hybridization,
};

/// Stable atom-table index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomId(usize);

impl AtomId {
    #[must_use]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for AtomId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// Stable bond-table index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondId(usize);

impl BondId {
    #[must_use]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }
    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for BondId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}
