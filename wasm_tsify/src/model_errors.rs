//! Exact detached-value errors used by molecular query projections.
use crate::errors::{enum_error, error_projection, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

enum_error!(StereoGroupError, StereoGroupErrorKind, ck::StereoGroupError, {
    DuplicateAtom => ck::StereoGroupError::DuplicateAtom,
    DuplicateBond => ck::StereoGroupError::DuplicateBond,
}, {});

enum_error!(TemplateAttachmentOrderError, TemplateAttachmentOrderErrorKind, ck::TemplateAttachmentOrderError, {
    Empty => ck::TemplateAttachmentOrderError::Empty,
    DuplicateTarget => ck::TemplateAttachmentOrderError::DuplicateTarget { .. },
    DuplicateLabel => ck::TemplateAttachmentOrderError::DuplicateLabel { .. },
    TargetOutOfRange => ck::TemplateAttachmentOrderError::TargetOutOfRange { .. },
    MappingTargetOutOfRange => ck::TemplateAttachmentOrderError::MappingTargetOutOfRange { .. },
    TargetRemoved => ck::TemplateAttachmentOrderError::TargetRemoved { .. },
}, {
    #[wasm_bindgen(getter, js_name = firstPosition)]
    pub fn first_position(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::TemplateAttachmentOrderError::DuplicateTarget { first_position, .. }
        | ck::TemplateAttachmentOrderError::DuplicateLabel { first_position, .. } => Some(*first_position), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = duplicatePosition)]
    pub fn duplicate_position(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::TemplateAttachmentOrderError::DuplicateTarget { duplicate_position, .. }
        | ck::TemplateAttachmentOrderError::DuplicateLabel { duplicate_position, .. } => Some(*duplicate_position), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn position(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::TemplateAttachmentOrderError::TargetOutOfRange { position, .. }
        | ck::TemplateAttachmentOrderError::MappingTargetOutOfRange { position, .. }
        | ck::TemplateAttachmentOrderError::TargetRemoved { position, .. } => Some(*position), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn target(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::TemplateAttachmentOrderError::DuplicateTarget { target, .. }
        | ck::TemplateAttachmentOrderError::TargetOutOfRange { target, .. }
        | ck::TemplateAttachmentOrderError::MappingTargetOutOfRange { target, .. }
        | ck::TemplateAttachmentOrderError::TargetRemoved { target, .. } => Some(target.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn label(&self) -> Option<String> { match self.inner.as_ref() {
        ck::TemplateAttachmentOrderError::DuplicateLabel { label, .. } => Some(label.clone()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = atomCount)]
    pub fn atom_count(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::TemplateAttachmentOrderError::TargetOutOfRange { atom_count, .. } => Some(*atom_count), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = mappingLen)]
    pub fn mapping_len(&self) -> Option<usize> { match self.inner.as_ref() {
        ck::TemplateAttachmentOrderError::MappingTargetOutOfRange { mapping_len, .. } => Some(*mapping_len), _ => None,
    } }
});

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "StereoGroupError | TemplateAttachmentOrderError")]
    #[derive(Clone)]
    pub type TopologyCause;
}

pub(crate) fn topology_cause(inner: &ck::TopologyValidationError) -> Option<TopologyCause> {
    match inner {
        ck::TopologyValidationError::StereoGroup(source) => {
            Some(StereoGroupError::from_rust(*source).unchecked_into())
        }
        ck::TopologyValidationError::TemplateAttachmentOrder { source, .. } => {
            Some(TemplateAttachmentOrderError::from_rust(source.clone()).unchecked_into())
        }
        ck::TopologyValidationError::AtomIdMismatch { .. }
        | ck::TopologyValidationError::BondIdMismatch { .. }
        | ck::TopologyValidationError::BondEndpointOutOfRange { .. }
        | ck::TopologyValidationError::SelfLoopBond { .. }
        | ck::TopologyValidationError::StereoAtomOutOfRange { .. }
        | ck::TopologyValidationError::StereoReferenceOutOfRange { .. }
        | ck::TopologyValidationError::StereoAtomsRequired { .. }
        | ck::TopologyValidationError::SubstanceGroupIdMismatch { .. }
        | ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { .. }
        | ck::TopologyValidationError::SubstanceGroupBondOutOfRange { .. }
        | ck::TopologyValidationError::SubstanceGroupParentOutOfRange { .. }
        | ck::TopologyValidationError::StereoGroupAtomOutOfRange { .. }
        | ck::TopologyValidationError::StereoGroupBondOutOfRange { .. }
        | ck::TopologyValidationError::AdjacencyMismatch => None,
    }
}

enum_error!(TopologyValidationError, TopologyValidationErrorKind, ck::TopologyValidationError, {
    StereoGroup => ck::TopologyValidationError::StereoGroup(_),
    AtomIdMismatch => ck::TopologyValidationError::AtomIdMismatch { .. },
    TemplateAttachmentOrder => ck::TopologyValidationError::TemplateAttachmentOrder { .. },
    BondIdMismatch => ck::TopologyValidationError::BondIdMismatch { .. },
    BondEndpointOutOfRange => ck::TopologyValidationError::BondEndpointOutOfRange { .. },
    SelfLoopBond => ck::TopologyValidationError::SelfLoopBond { .. },
    StereoAtomOutOfRange => ck::TopologyValidationError::StereoAtomOutOfRange { .. },
    StereoReferenceOutOfRange => ck::TopologyValidationError::StereoReferenceOutOfRange { .. },
    StereoAtomsRequired => ck::TopologyValidationError::StereoAtomsRequired { .. },
    SubstanceGroupIdMismatch => ck::TopologyValidationError::SubstanceGroupIdMismatch { .. },
    SubstanceGroupAtomOutOfRange => ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { .. },
    SubstanceGroupBondOutOfRange => ck::TopologyValidationError::SubstanceGroupBondOutOfRange { .. },
    SubstanceGroupParentOutOfRange => ck::TopologyValidationError::SubstanceGroupParentOutOfRange { .. },
    StereoGroupAtomOutOfRange => ck::TopologyValidationError::StereoGroupAtomOutOfRange { .. },
    StereoGroupBondOutOfRange => ck::TopologyValidationError::StereoGroupBondOutOfRange { .. },
    AdjacencyMismatch => ck::TopologyValidationError::AdjacencyMismatch,
}, {
    #[wasm_bindgen(getter)]
    pub fn position(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::AtomIdMismatch { position, .. } | ck::TopologyValidationError::BondIdMismatch { position, .. }
        | ck::TopologyValidationError::SubstanceGroupIdMismatch { position, .. } => Some(position), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::AtomIdMismatch { id, .. } => Some(id.index()),
        ck::TopologyValidationError::BondIdMismatch { id, .. } => Some(id.index()),
        ck::TopologyValidationError::SubstanceGroupIdMismatch { id, .. } => Some(id.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::TemplateAttachmentOrder { atom, .. } | ck::TopologyValidationError::BondEndpointOutOfRange { atom, .. }
        | ck::TopologyValidationError::SelfLoopBond { atom, .. } | ck::TopologyValidationError::StereoReferenceOutOfRange { atom, .. }
        | ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { atom, .. } | ck::TopologyValidationError::StereoGroupAtomOutOfRange { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::BondEndpointOutOfRange { bond, .. } | ck::TopologyValidationError::SelfLoopBond { bond, .. }
        | ck::TopologyValidationError::StereoAtomOutOfRange { bond, .. } | ck::TopologyValidationError::StereoReferenceOutOfRange { bond, .. }
        | ck::TopologyValidationError::StereoAtomsRequired { bond, .. } | ck::TopologyValidationError::SubstanceGroupBondOutOfRange { bond, .. }
        | ck::TopologyValidationError::StereoGroupBondOutOfRange { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = atomCount)]
    pub fn atom_count(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::BondEndpointOutOfRange { atom_count, .. } | ck::TopologyValidationError::StereoAtomOutOfRange { atom_count, .. }
        | ck::TopologyValidationError::StereoReferenceOutOfRange { atom_count, .. } | ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { atom_count, .. }
        | ck::TopologyValidationError::StereoGroupAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = bondCount)]
    pub fn bond_count(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::SubstanceGroupBondOutOfRange { bond_count, .. }
        | ck::TopologyValidationError::StereoGroupBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn endpoint(&self) -> Option<String> { match *self.inner {
        ck::TopologyValidationError::BondEndpointOutOfRange { endpoint, .. } => Some(endpoint.into()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn begin(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::StereoAtomOutOfRange { begin, .. } => Some(begin.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::StereoAtomOutOfRange { end, .. } => Some(end.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn stereo(&self) -> Option<crate::enums::BondStereo> { match *self.inner {
        ck::TopologyValidationError::StereoAtomsRequired { stereo, .. } => Some(stereo.into()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn sgroup(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::SubstanceGroupAtomOutOfRange { sgroup, .. }
        | ck::TopologyValidationError::SubstanceGroupBondOutOfRange { sgroup, .. }
        | ck::TopologyValidationError::SubstanceGroupParentOutOfRange { sgroup, .. } => Some(sgroup.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn parent(&self) -> Option<usize> { match *self.inner {
        ck::TopologyValidationError::SubstanceGroupParentOutOfRange { parent, .. } => Some(parent.index()), _ => None,
    } }
}, source = TopologyCause, cause = topology_cause);

enum_error!(AdjacencyError, AdjacencyErrorKind, ck::AdjacencyError, {
    BondAtomOutOfRange => ck::AdjacencyError::BondAtomOutOfRange { .. },
    DuplicateBondId => ck::AdjacencyError::DuplicateBondId { .. },
    DuplicateEdge => ck::AdjacencyError::DuplicateEdge { .. },
}, {
    #[wasm_bindgen(getter)]
    pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::BondAtomOutOfRange { bond, .. } | ck::AdjacencyError::DuplicateBondId { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::BondAtomOutOfRange { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn endpoint(&self) -> Option<String> { match *self.inner {
        ck::AdjacencyError::BondAtomOutOfRange { endpoint, .. } => Some(endpoint.into()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = atomCount)]
    pub fn atom_count(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::BondAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = firstPosition)]
    pub fn first_position(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::DuplicateBondId { first_position, .. } => Some(first_position), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = secondPosition)]
    pub fn second_position(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::DuplicateBondId { second_position, .. } => Some(second_position), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = firstBond)]
    pub fn first_bond(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::DuplicateEdge { first_bond, .. } => Some(first_bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = secondBond)]
    pub fn second_bond(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::DuplicateEdge { second_bond, .. } => Some(second_bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn begin(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::DuplicateEdge { begin, .. } => Some(begin.index()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> Option<usize> { match *self.inner {
        ck::AdjacencyError::DuplicateEdge { end, .. } => Some(end.index()), _ => None,
    } }
});

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum PropertyValueKind {
    String,
    Int,
    UInt,
    IntVector,
    StringVector,
    Double,
    Bool,
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "PropertyValueError")]
    pub type PropertyValueErrorValue;
}
impl From<ck::PropertyValueKind> for PropertyValueKind {
    fn from(value: ck::PropertyValueKind) -> Self {
        match value {
            ck::PropertyValueKind::String => Self::String,
            ck::PropertyValueKind::Int => Self::Int,
            ck::PropertyValueKind::UInt => Self::UInt,
            ck::PropertyValueKind::IntVector => Self::IntVector,
            ck::PropertyValueKind::StringVector => Self::StringVector,
            ck::PropertyValueKind::Double => Self::Double,
            ck::PropertyValueKind::Bool => Self::Bool,
        }
    }
}
error_projection!(PropertyValueError, ck::PropertyValueError, {
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> PropertyValueKind {
        self.inner.expected().into()
    }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> PropertyValueKind {
        self.inner.actual().into()
    }
});
enum_error!(MoleculePropertyError, MoleculePropertyErrorKind, ck::MoleculePropertyError, {
    EmptyKey => ck::MoleculePropertyError::EmptyKey,
    ComputedListKind => ck::MoleculePropertyError::ComputedListKind(_),
}, {
    #[wasm_bindgen(getter, js_name = computedListKind)]
    pub fn computed_list_kind(&self) -> Option<PropertyValueErrorValue> { match self.inner.as_ref() {
        ck::MoleculePropertyError::ComputedListKind(value) => Some(PropertyValueError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
});

enum_error!(AtomPropertyError, AtomPropertyErrorKind, ck::AtomPropertyError, {
    EmptyKey => ck::AtomPropertyError::EmptyKey,
    ComputedListKind => ck::AtomPropertyError::ComputedListKind(_),
}, {
    #[wasm_bindgen(getter, js_name = computedListKind)]
    pub fn computed_list_kind(&self) -> Option<PropertyValueErrorValue> { match self.inner.as_ref() {
        ck::AtomPropertyError::ComputedListKind(value) => Some(PropertyValueError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
});
enum_error!(BondValueError, BondValueErrorKind, ck::BondValueError, {
    StereoAtomsRequired => ck::BondValueError::StereoAtomsRequired,
    EmptyPropertyKey => ck::BondValueError::EmptyPropertyKey,
    ComputedListKind => ck::BondValueError::ComputedListKind(_),
}, {
    #[wasm_bindgen(getter, js_name = computedListKind)]
    pub fn computed_list_kind(&self) -> Option<PropertyValueErrorValue> { match self.inner.as_ref() {
        ck::BondValueError::ComputedListKind(value) => Some(PropertyValueError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
});

enum_error!(MappingValidationError, MappingValidationErrorKind, ck::MappingValidationError, {
    Length => ck::MappingValidationError::Length { .. },
    OutOfRange => ck::MappingValidationError::OutOfRange { .. },
    InverseMismatch => ck::MappingValidationError::InverseMismatch { .. },
}, {
    #[wasm_bindgen(getter)]
    pub fn entity(&self) -> String { match *self.inner {
        ck::MappingValidationError::Length { entity, .. } | ck::MappingValidationError::OutOfRange { entity, .. }
        | ck::MappingValidationError::InverseMismatch { entity, .. } => entity.into(),
    } }
    #[wasm_bindgen(getter)]
    pub fn direction(&self) -> String { match *self.inner {
        ck::MappingValidationError::Length { direction, .. } | ck::MappingValidationError::OutOfRange { direction, .. }
        | ck::MappingValidationError::InverseMismatch { direction, .. } => direction.into(),
    } }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> Option<usize> { match *self.inner { ck::MappingValidationError::Length { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> Option<usize> { match *self.inner { ck::MappingValidationError::Length { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn row(&self) -> Option<usize> { match *self.inner {
        ck::MappingValidationError::OutOfRange { row, .. } | ck::MappingValidationError::InverseMismatch { row, .. } => Some(row), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn mapped(&self) -> Option<usize> { match *self.inner {
        ck::MappingValidationError::OutOfRange { mapped, .. } | ck::MappingValidationError::InverseMismatch { mapped, .. } => Some(mapped), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = targetCount)]
    pub fn target_count(&self) -> Option<usize> { match *self.inner { ck::MappingValidationError::OutOfRange { target_count, .. } => Some(target_count), _ => None } }
});

enum_error!(CoordinateValidationError, CoordinateValidationErrorKind, ck::CoordinateValidationError, {
    AtomIndexOverflow => ck::CoordinateValidationError::AtomIndexOverflow { .. },
    Missing3DConformer => ck::CoordinateValidationError::Missing3DConformer { .. },
    ConformerIdOverflow => ck::CoordinateValidationError::ConformerIdOverflow { .. },
    RowCount => ck::CoordinateValidationError::RowCount { .. },
    MissingSourceConformerOrder => ck::CoordinateValidationError::MissingSourceConformerOrder,
    SourceConformerOrder => ck::CoordinateValidationError::SourceConformerOrder { .. },
    DuplicateConformerId => ck::CoordinateValidationError::DuplicateConformerId { .. },
    NonFiniteCoordinate => ck::CoordinateValidationError::NonFiniteCoordinate { .. },
}, {
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::CoordinateValidationError::AtomIndexOverflow { atom } | ck::CoordinateValidationError::NonFiniteCoordinate { atom, .. } => Some(atom), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn id(&self) -> Option<usize> { match *self.inner {
        ck::CoordinateValidationError::Missing3DConformer { id } | ck::CoordinateValidationError::DuplicateConformerId { id, .. } => Some(id), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = maxId)]
    pub fn max_id(&self) -> Option<usize> { match *self.inner { ck::CoordinateValidationError::ConformerIdOverflow { max_id } => Some(max_id), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn dimension(&self) -> Option<String> { match *self.inner {
        ck::CoordinateValidationError::RowCount { dimension, .. } | ck::CoordinateValidationError::DuplicateConformerId { dimension, .. }
        | ck::CoordinateValidationError::NonFiniteCoordinate { dimension, .. } => Some(dimension.into()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn conformer(&self) -> Option<usize> { match *self.inner {
        ck::CoordinateValidationError::RowCount { conformer, .. } | ck::CoordinateValidationError::NonFiniteCoordinate { conformer, .. } => Some(conformer), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn rows(&self) -> Option<usize> { match *self.inner { ck::CoordinateValidationError::RowCount { rows, .. } => Some(rows), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)]
    pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::CoordinateValidationError::RowCount { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = twoD)]
    pub fn two_d(&self) -> Option<usize> { match *self.inner { ck::CoordinateValidationError::SourceConformerOrder { two_d, .. } => Some(two_d), _ => None } }
    #[wasm_bindgen(getter, js_name = threeD)]
    pub fn three_d(&self) -> Option<usize> { match *self.inner { ck::CoordinateValidationError::SourceConformerOrder { three_d, .. } => Some(three_d), _ => None } }
    #[wasm_bindgen(getter, js_name = expectedTwoD)]
    pub fn expected_two_d(&self) -> Option<usize> { match *self.inner { ck::CoordinateValidationError::SourceConformerOrder { expected_two_d, .. } => Some(expected_two_d), _ => None } }
    #[wasm_bindgen(getter, js_name = expectedThreeD)]
    pub fn expected_three_d(&self) -> Option<usize> { match *self.inner { ck::CoordinateValidationError::SourceConformerOrder { expected_three_d, .. } => Some(expected_three_d), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn axis(&self) -> Option<String> { match *self.inner { ck::CoordinateValidationError::NonFiniteCoordinate { axis, .. } => Some(axis.into()), _ => None } }
});

enum_error!(CipDescriptorError, CipDescriptorErrorKind, ck::CipDescriptorError, {
    InvalidNeighborOrder => ck::CipDescriptorError::InvalidNeighborOrder { .. },
    Property => ck::CipDescriptorError::Property(_),
    InvalidStoredDescriptor => ck::CipDescriptorError::InvalidStoredDescriptor { .. },
}, {
    #[wasm_bindgen(getter)]
    pub fn value(&self) -> Result<Option<String>, JsValue> { match self.inner.as_ref() {
        ck::CipDescriptorError::InvalidNeighborOrder { value, .. } | ck::CipDescriptorError::InvalidStoredDescriptor { value } =>
            std::str::from_utf8(value.as_bytes()).map(|value| Some(value.to_owned())).map_err(|error| js_sys::TypeError::new(&error.to_string()).into()),
        ck::CipDescriptorError::Property(_) => Ok(None),
    } }
    #[wasm_bindgen(getter)]
    pub fn detail(&self) -> Option<String> { match self.inner.as_ref() {
        ck::CipDescriptorError::InvalidNeighborOrder { detail, .. } => Some(detail.clone()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn property(&self) -> Option<PropertyValueErrorValue> { match self.inner.as_ref() {
        ck::CipDescriptorError::Property(value) => Some(PropertyValueError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
});
