//! Detached edit errors: typed payloads are not invented source links.
use crate::{
    errors::{enum_error, exception},
    model_errors::*,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum BondEndPointsParseErrorKind {
    InvalidArgument,
    OutOfRange,
    UnrepresentableEndpoint,
}
impl From<ck::BondEndPointsParseErrorKind> for BondEndPointsParseErrorKind {
    fn from(value: ck::BondEndPointsParseErrorKind) -> Self {
        match value {
            ck::BondEndPointsParseErrorKind::InvalidArgument => Self::InvalidArgument,
            ck::BondEndPointsParseErrorKind::OutOfRange => Self::OutOfRange,
            ck::BondEndPointsParseErrorKind::UnrepresentableEndpoint => {
                Self::UnrepresentableEndpoint
            }
        }
    }
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | CoordinateValidationError | MoleculePropertyError | AtomPropertyError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type EditCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type TopologyPayload;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type BondPayload;
    #[wasm_bindgen(typescript_type = "MappingValidationError")]
    pub type MappingPayload;
}
pub(crate) fn cause(inner: &ck::TopologyEditError) -> Option<EditCause> {
    match inner {
        ck::TopologyEditError::StereoGroup(value) => {
            Some(StereoGroupError::from_rust(*value).unchecked_into())
        }
        ck::TopologyEditError::Coordinate(value) => {
            Some(CoordinateValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::TopologyEditError::MoleculeProperty(value) => {
            Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::TopologyEditError::AtomProperty(value) => {
            Some(AtomPropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::TopologyEditError::TemplateAttachmentRemap { source, .. } => {
            Some(TemplateAttachmentOrderError::from_rust(source.clone()).unchecked_into())
        }
        ck::TopologyEditError::SourceGroupParentCycle { .. }
        | ck::TopologyEditError::InvalidSource(_)
        | ck::TopologyEditError::AtomOutOfRange { .. }
        | ck::TopologyEditError::BondOutOfRange { .. }
        | ck::TopologyEditError::BondEndPointsParse { .. }
        | ck::TopologyEditError::IncompleteSourceBatchMasks
        | ck::TopologyEditError::NullBondAtom
        | ck::TopologyEditError::DuplicateBond { .. }
        | ck::TopologyEditError::InvalidBond(_)
        | ck::TopologyEditError::PermutationLength { .. }
        | ck::TopologyEditError::PermutationAtomOutOfRange { .. }
        | ck::TopologyEditError::PermutationDuplicateAtom { .. }
        | ck::TopologyEditError::InvalidMapping(_)
        | ck::TopologyEditError::InvalidResult(_) => None,
    }
}
enum_error!(TopologyEditError, TopologyEditErrorKind, ck::TopologyEditError, {
    StereoGroup => ck::TopologyEditError::StereoGroup(_),
    Coordinate => ck::TopologyEditError::Coordinate(_),
    MoleculeProperty => ck::TopologyEditError::MoleculeProperty(_),
    AtomProperty => ck::TopologyEditError::AtomProperty(_),
    SourceGroupParentCycle => ck::TopologyEditError::SourceGroupParentCycle { .. },
    InvalidSource => ck::TopologyEditError::InvalidSource(_),
    AtomOutOfRange => ck::TopologyEditError::AtomOutOfRange { .. },
    BondOutOfRange => ck::TopologyEditError::BondOutOfRange { .. },
    BondEndPointsParse => ck::TopologyEditError::BondEndPointsParse { .. },
    IncompleteSourceBatchMasks => ck::TopologyEditError::IncompleteSourceBatchMasks,
    NullBondAtom => ck::TopologyEditError::NullBondAtom,
    DuplicateBond => ck::TopologyEditError::DuplicateBond { .. },
    InvalidBond => ck::TopologyEditError::InvalidBond(_),
    PermutationLength => ck::TopologyEditError::PermutationLength { .. },
    PermutationAtomOutOfRange => ck::TopologyEditError::PermutationAtomOutOfRange { .. },
    PermutationDuplicateAtom => ck::TopologyEditError::PermutationDuplicateAtom { .. },
    InvalidMapping => ck::TopologyEditError::InvalidMapping(_),
    TemplateAttachmentRemap => ck::TopologyEditError::TemplateAttachmentRemap { .. },
    InvalidResult => ck::TopologyEditError::InvalidResult(_),
}, {
    #[wasm_bindgen(getter)]
    pub fn group(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::SourceGroupParentCycle { group } => Some(group.index()), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::TopologyEditError::AtomOutOfRange { atom, .. } | ck::TopologyEditError::PermutationAtomOutOfRange { atom, .. }
        | ck::TopologyEditError::PermutationDuplicateAtom { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = atomCount)]
    pub fn atom_count(&self) -> Option<usize> { match *self.inner {
        ck::TopologyEditError::AtomOutOfRange { atom_count, .. } | ck::TopologyEditError::PermutationAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::TopologyEditError::BondOutOfRange { bond, .. } | ck::TopologyEditError::BondEndPointsParse { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = bondCount)]
    pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::BondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter, js_name = tokenIndex)]
    pub fn token_index(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::BondEndPointsParse { token_index, .. } => Some(token_index), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn token(&self) -> Result<Option<String>, JsValue> { match self.inner.as_ref() {
        ck::TopologyEditError::BondEndPointsParse { token, .. } => std::str::from_utf8(token.as_bytes()).map(|value| Some(value.to_owned())).map_err(|error| js_sys::TypeError::new(&error.to_string()).into()),
        _ => Ok(None),
    } }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> Option<BondEndPointsParseErrorKind> { match *self.inner { ck::TopologyEditError::BondEndPointsParse { kind, .. } => Some(kind.into()), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn begin(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::DuplicateBond { begin, .. } => Some(begin.index()), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn end(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::DuplicateBond { end, .. } => Some(end.index()), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::PermutationLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::PermutationLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)]
    pub fn position(&self) -> Option<usize> { match *self.inner {
        ck::TopologyEditError::PermutationAtomOutOfRange { position, .. } | ck::TopologyEditError::PermutationDuplicateAtom { position, .. } => Some(position), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn carrier(&self) -> Option<usize> { match *self.inner { ck::TopologyEditError::TemplateAttachmentRemap { carrier, .. } => Some(carrier.index()), _ => None } }
    #[wasm_bindgen(getter, js_name = invalidSource)]
    pub fn invalid_source(&self) -> Option<TopologyPayload> { match self.inner.as_ref() {
        ck::TopologyEditError::InvalidSource(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = invalidResult)]
    pub fn invalid_result(&self) -> Option<TopologyPayload> { match self.inner.as_ref() {
        ck::TopologyEditError::InvalidResult(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = invalidBond)]
    pub fn invalid_bond(&self) -> Option<BondPayload> { match self.inner.as_ref() {
        ck::TopologyEditError::InvalidBond(value) => Some(BondValueError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = invalidMapping)]
    pub fn invalid_mapping(&self) -> Option<MappingPayload> { match self.inner.as_ref() {
        ck::TopologyEditError::InvalidMapping(value) => Some(MappingValidationError::from_rust(value.clone()).unchecked_into()), _ => None,
    } }
}, source = EditCause, cause = cause);
