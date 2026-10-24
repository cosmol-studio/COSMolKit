//! Ring errors preserve native transparent-source versus payload distinctions.
use crate::{
    errors::{enum_error, exception},
    model_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

enum_error!(RingDecomposerError, RingDecomposerErrorKind, ck::RingDecomposerError, {
    NodeOutOfRange => ck::RingDecomposerError::NodeOutOfRange { .. },
    SelfLoop => ck::RingDecomposerError::SelfLoop { .. },
    DuplicateEdge => ck::RingDecomposerError::DuplicateEdge { .. },
    EmptyGraph => ck::RingDecomposerError::EmptyGraph,
    EdgeNotFound => ck::RingDecomposerError::EdgeNotFound { .. },
    RelevantCycleCountOverflow => ck::RingDecomposerError::RelevantCycleCountOverflow,
    InvalidShortestPathDag => ck::RingDecomposerError::InvalidShortestPathDag { .. },
}, {
    #[wasm_bindgen(getter)]
    pub fn node(&self) -> Option<usize> { match *self.inner {
        ck::RingDecomposerError::NodeOutOfRange { node, .. } | ck::RingDecomposerError::SelfLoop { node } => Some(node), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = nodeCount)]
    pub fn node_count(&self) -> Option<usize> { match *self.inner {
        ck::RingDecomposerError::NodeOutOfRange { node_count, .. } => Some(node_count), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn from(&self) -> Option<usize> { match *self.inner {
        ck::RingDecomposerError::DuplicateEdge { from, .. } | ck::RingDecomposerError::EdgeNotFound { from, .. } => Some(from), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn to(&self) -> Option<usize> { match *self.inner {
        ck::RingDecomposerError::DuplicateEdge { to, .. } | ck::RingDecomposerError::EdgeNotFound { to, .. } => Some(to), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn root(&self) -> Option<usize> { match *self.inner {
        ck::RingDecomposerError::InvalidShortestPathDag { root, .. } => Some(root), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn target(&self) -> Option<usize> { match *self.inner {
        ck::RingDecomposerError::InvalidShortestPathDag { target, .. } => Some(target), _ => None,
    } }
});

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "MoleculePropertyError | StereoGroupError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type RingCause;
    #[wasm_bindgen(typescript_type = "AdjacencyError")]
    pub type AdjacencyPayload;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type TopologyPayload;
    #[wasm_bindgen(typescript_type = "RingDecomposerError")]
    pub type DecomposerPayload;
}
pub(crate) fn cause(inner: &ck::RingFindingError) -> Option<RingCause> {
    match inner {
        ck::RingFindingError::MoleculeProperty(value) => {
            Some(model_errors::MoleculePropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::RingFindingError::InvalidTopology(value) => {
            model_errors::topology_cause(value).map(JsCast::unchecked_into)
        }
        ck::RingFindingError::SourceRingFindType { .. }
        | ck::RingFindingError::Value { .. }
        | ck::RingFindingError::ExpectedBondNotFound { .. }
        | ck::RingFindingError::RingAtomOutOfRange { .. }
        | ck::RingFindingError::RingBondOutOfRange { .. }
        | ck::RingFindingError::RingFamilyEdgeMappingMissing { .. }
        | ck::RingFindingError::RelevantCycleCountOutOfRange
        | ck::RingFindingError::UnsupportedBranch { .. }
        | ck::RingFindingError::Adjacency(_)
        | ck::RingFindingError::RingDecomposer(_) => None,
    }
}
enum_error!(RingFindingError, RingFindingErrorKind, ck::RingFindingError, {
    SourceRingFindType => ck::RingFindingError::SourceRingFindType { .. },
    MoleculeProperty => ck::RingFindingError::MoleculeProperty(_),
    Value => ck::RingFindingError::Value { .. },
    ExpectedBondNotFound => ck::RingFindingError::ExpectedBondNotFound { .. },
    RingAtomOutOfRange => ck::RingFindingError::RingAtomOutOfRange { .. },
    RingBondOutOfRange => ck::RingFindingError::RingBondOutOfRange { .. },
    RingFamilyEdgeMappingMissing => ck::RingFindingError::RingFamilyEdgeMappingMissing { .. },
    RelevantCycleCountOutOfRange => ck::RingFindingError::RelevantCycleCountOutOfRange,
    UnsupportedBranch => ck::RingFindingError::UnsupportedBranch { .. },
    Adjacency => ck::RingFindingError::Adjacency(_),
    InvalidTopology => ck::RingFindingError::InvalidTopology(_),
    RingDecomposer => ck::RingFindingError::RingDecomposer(_),
}, {
    #[wasm_bindgen(getter)]
    pub fn tag(&self) -> Option<u32> { match *self.inner { ck::RingFindingError::SourceRingFindType { tag } => Some(tag), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn begin(&self) -> Option<usize> { match *self.inner { ck::RingFindingError::ExpectedBondNotFound { begin, .. } => Some(begin.index()), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> Option<usize> { match *self.inner { ck::RingFindingError::ExpectedBondNotFound { end, .. } => Some(end.index()), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> Option<usize> { match *self.inner { ck::RingFindingError::RingAtomOutOfRange { atom, .. } => Some(atom), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)]
    pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::RingFindingError::RingAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn bond(&self) -> Option<usize> { match *self.inner { ck::RingFindingError::RingBondOutOfRange { bond, .. } => Some(bond), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)]
    pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::RingFindingError::RingBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn edge(&self) -> Option<usize> { match *self.inner { ck::RingFindingError::RingFamilyEdgeMappingMissing { edge } => Some(edge), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn reason(&self) -> Option<String> { match *self.inner { ck::RingFindingError::UnsupportedBranch { reason } => Some(reason.into()), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn adjacency(&self) -> Option<AdjacencyPayload> { match self.inner.as_ref() {
        ck::RingFindingError::Adjacency(value) => Some(model_errors::AdjacencyError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn topology(&self) -> Option<TopologyPayload> { match self.inner.as_ref() {
        ck::RingFindingError::InvalidTopology(value) => Some(model_errors::TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = ringDecomposer)]
    pub fn ring_decomposer(&self) -> Option<DecomposerPayload> { match self.inner.as_ref() {
        ck::RingFindingError::RingDecomposer(value) => Some(RingDecomposerError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
}, source = RingCause, cause = cause);
