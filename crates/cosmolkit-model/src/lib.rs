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
mod mapping;
mod properties;
mod query;
mod sgroup;
mod topology;

pub use adjacency::{AdjacencyError, AdjacencyList, NeighborRef};
pub use atom::{
    Atom, AtomId, AtomPdbResidueInfo, AtomPropertyError, AtomSpec, TemplateAttachment,
    TemplateAttachmentOrder, TemplateAttachmentOrderError,
};
pub use bond::{Bond, BondId, BondSpec, BondValueError};
pub use cip::{CipDescriptor, CipDescriptorError};
pub use coordinates::{
    Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateValidationError,
};
pub use mapping::{AtomMapping, BondMapping, MappingValidationError, TopologyMapping};
pub use properties::{
    MoleculeProperties, MoleculePropertyError, SdfPropertyList, SdfPropertyListTarget,
};
pub use query::{
    AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, BondQueryPredicate,
    QueryAtom, QueryBond, QueryGraph, QueryGraphError, QueryNode, QueryStateError, QueryStateRef,
    RecursiveStructureQuery, remap_query_rows,
};
pub use sgroup::{
    SGroupAttachPoint, SGroupBondRole, SGroupBracket, SGroupBracketStyle, SGroupCState,
    SGroupConnection, SGroupData, SGroupDisplay, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupId, SubstanceGroupKind,
};
pub use topology::{TopologyBatchEdit, TopologyBlock, TopologyEditError, TopologyValidationError};

pub use cosmolkit_types::{
    BondDirection, BondOrder, BondStereo, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element,
    ElementInfo, ElementParseError, Hybridization,
};
