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
mod property_text;
mod property_value;
mod query;
mod sgroup;
mod source_atom_facts;
mod stereo;
mod topology;
mod valence;

pub use adjacency::{AdjacencyError, AdjacencyList, NeighborRef};
pub use atom::{
    Atom, AtomId, AtomPdbResidueInfo, AtomPropertyError, AtomSpec, TemplateAttachment,
    TemplateAttachmentOrder, TemplateAttachmentOrderError, ordered_atom_properties,
    replace_atom_template_attachment_order,
};
pub use bond::{Bond, BondId, BondSpec, BondValueError, ordered_bond_properties};
pub use cip::{CipDescriptor, CipDescriptorError};
pub use coordinates::{
    Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, CoordinateSourceConformer,
    CoordinateValidationError, first_non_finite_coordinate,
};
pub use mapping::{AtomMapping, BondMapping, MappingValidationError, TopologyMapping};
pub use properties::{
    MoleculeProperties, MoleculePropertyError, SdfPropertyList, SdfPropertyListTarget,
};
#[doc(hidden)]
pub use property_text::PropertyText;
#[doc(hidden)]
pub use property_value::{MissingPropertyError, PropertyStore, PropertyStoreError};
pub use property_value::{PropertyValue, PropertyValueError, PropertyValueKind};
pub use query::{
    AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, BondQueryPredicate,
    QueryAtom, QueryAtomConversionError, QueryAtomIdentity, QueryBond, QueryGraph, QueryGraphError,
    QueryNode, QueryStateError, QueryStateRef, RecursiveStructureQuery,
    ordered_query_atom_properties, query_substance_groups, remap_query_rows,
    remap_query_rows_with_appended, replace_query_stereo_groups, replace_query_substance_groups,
};
pub use sgroup::{
    SGroupAttachPoint, SGroupBondRole, SGroupBracket, SGroupBracketStyle, SGroupCState,
    SGroupConnection, SGroupData, SGroupDisplay, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupId, SubstanceGroupKind, set_stereo_group_write_id, stereo_group_write_id,
};
pub use topology::{
    BondEndPointsParseErrorKind, TopologyBatchEdit, TopologyBlock, TopologyEditError,
    TopologyValidationError,
};

pub use cosmolkit_types::{
    BondDirection, BondOrder, BondStereo, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element,
    ElementInfo, ElementParseError, Hybridization,
};

pub use valence::{AtomMetadata, ValenceError, ValencePhase};

#[doc(hidden)]
pub use source_atom_facts::SourceAtomValenceFacts;

pub use stereo::{LigandRef, TetrahedralStereo};

#[doc(hidden)]
pub use sgroup::merge_absolute_stereo_groups;

#[doc(hidden)]
pub use sgroup::insert_stereo_groups;

#[doc(hidden)]
pub use topology::{SourceBondBatchMasks, SourceBondNeighbors, add_source_bond_order};

// Narrow source state/value boundary for detached owning algorithms only.
#[doc(hidden)]
pub use topology::{SourceBatchCommitState, commit_batch_edit_source};

#[doc(hidden)]
pub use topology::{replace_source_bond, source_bond_between_atoms};

#[doc(hidden)]
pub use topology::add_source_bond_value;

#[doc(hidden)]
pub use coordinates::source_set_atom_position;

mod source_ring_info;
#[doc(hidden)]
pub use source_ring_info::SourceRingInfo;
