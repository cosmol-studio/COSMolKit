//! Canonical ranking errors retain source-defined transparent edges.
use crate::{
    errors::{enum_error, exception},
    model_errors::{PropertyValueKind, StereoGroupError, TemplateAttachmentOrderError},
    ring_errors, valence_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | TemplateAttachmentOrderError | MoleculePropertyError | TopologyValidationError"
    )]
    #[derive(Clone)]
    pub type CanonicalCause;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type CanonicalRingsPayload;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type CanonicalValencePayload;
}
pub(crate) fn cause(value: &ck::CanonicalRankError) -> Option<CanonicalCause> {
    match value {
        ck::CanonicalRankError::StereoGroup(value) => {
            Some(StereoGroupError::from_rust(*value).unchecked_into())
        }
        ck::CanonicalRankError::TemplateAttachmentRemap { source, .. } => {
            Some(TemplateAttachmentOrderError::from_rust(source.clone()).unchecked_into())
        }
        ck::CanonicalRankError::RingFinding(value) => {
            ring_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CanonicalRankError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CanonicalRankError::HanoiScratchTooSmall
        | ck::CanonicalRankError::UnsignedRankOverflow { .. }
        | ck::CanonicalRankError::InvalidPropertyKind { .. }
        | ck::CanonicalRankError::ComponentAtomOutOfRange { .. }
        | ck::CanonicalRankError::DuplicateComponentAtom { .. }
        | ck::CanonicalRankError::AtomMaskLength { .. }
        | ck::CanonicalRankError::BondMaskLength { .. }
        | ck::CanonicalRankError::AtomSymbolLength { .. }
        | ck::CanonicalRankError::BondSymbolLength { .. }
        | ck::CanonicalRankError::PreparedValenceLength { .. }
        | ck::CanonicalRankError::PreparedValenceInvalid { .. }
        | ck::CanonicalRankError::PreparedRingLength { .. }
        | ck::CanonicalRankError::ProtocolDebt { .. } => None,
    }
}
enum_error!(CanonicalRankError, CanonicalRankErrorKind, ck::CanonicalRankError, {
    HanoiScratchTooSmall => ck::CanonicalRankError::HanoiScratchTooSmall, StereoGroup => ck::CanonicalRankError::StereoGroup(_),
    UnsignedRankOverflow => ck::CanonicalRankError::UnsignedRankOverflow { .. }, InvalidPropertyKind => ck::CanonicalRankError::InvalidPropertyKind { .. },
    ComponentAtomOutOfRange => ck::CanonicalRankError::ComponentAtomOutOfRange { .. }, DuplicateComponentAtom => ck::CanonicalRankError::DuplicateComponentAtom { .. },
    AtomMaskLength => ck::CanonicalRankError::AtomMaskLength { .. }, BondMaskLength => ck::CanonicalRankError::BondMaskLength { .. },
    AtomSymbolLength => ck::CanonicalRankError::AtomSymbolLength { .. }, BondSymbolLength => ck::CanonicalRankError::BondSymbolLength { .. },
    PreparedValenceLength => ck::CanonicalRankError::PreparedValenceLength { .. }, PreparedValenceInvalid => ck::CanonicalRankError::PreparedValenceInvalid { .. },
    PreparedRingLength => ck::CanonicalRankError::PreparedRingLength { .. }, ProtocolDebt => ck::CanonicalRankError::ProtocolDebt { .. },
    TemplateAttachmentRemap => ck::CanonicalRankError::TemplateAttachmentRemap { .. }, RingFinding => ck::CanonicalRankError::RingFinding(_), Valence => ck::CanonicalRankError::Valence(_),
}, {
    #[wasm_bindgen(getter, js_name = atomIndex)] pub fn atom_index(&self) -> Option<usize> { match *self.inner {
        ck::CanonicalRankError::UnsignedRankOverflow { atom_index, .. } | ck::CanonicalRankError::InvalidPropertyKind { atom_index, .. }
        | ck::CanonicalRankError::ComponentAtomOutOfRange { atom_index, .. } | ck::CanonicalRankError::DuplicateComponentAtom { atom_index } | ck::CanonicalRankError::PreparedValenceInvalid { atom_index } => Some(atom_index), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn property(&self) -> Option<String> { match *self.inner { ck::CanonicalRankError::UnsignedRankOverflow { property, .. } | ck::CanonicalRankError::InvalidPropertyKind { property, .. } => Some(property.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<u32> { match *self.inner { ck::CanonicalRankError::UnsignedRankOverflow { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> Option<PropertyValueKind> { match *self.inner { ck::CanonicalRankError::InvalidPropertyKind { kind, .. } => Some(kind.into()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::ComponentAtomOutOfRange { atom_count, .. } | ck::CanonicalRankError::PreparedValenceLength { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::AtomMaskLength { expected, .. } | ck::CanonicalRankError::BondMaskLength { expected, .. } | ck::CanonicalRankError::AtomSymbolLength { expected, .. } | ck::CanonicalRankError::BondSymbolLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::AtomMaskLength { actual, .. } | ck::CanonicalRankError::BondMaskLength { actual, .. } | ck::CanonicalRankError::AtomSymbolLength { actual, .. } | ck::CanonicalRankError::BondSymbolLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter, js_name = explicitLen)] pub fn explicit_len(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::PreparedValenceLength { explicit_len, .. } => Some(explicit_len), _ => None } }
    #[wasm_bindgen(getter, js_name = implicitLen)] pub fn implicit_len(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::PreparedValenceLength { implicit_len, .. } => Some(implicit_len), _ => None } }
    #[wasm_bindgen(getter, js_name = expectedAtoms)] pub fn expected_atoms(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::PreparedRingLength { expected_atoms, .. } => Some(expected_atoms), _ => None } }
    #[wasm_bindgen(getter, js_name = actualAtoms)] pub fn actual_atoms(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::PreparedRingLength { actual_atoms, .. } => Some(actual_atoms), _ => None } }
    #[wasm_bindgen(getter, js_name = expectedBonds)] pub fn expected_bonds(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::PreparedRingLength { expected_bonds, .. } => Some(expected_bonds), _ => None } }
    #[wasm_bindgen(getter, js_name = actualBonds)] pub fn actual_bonds(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::PreparedRingLength { actual_bonds, .. } => Some(actual_bonds), _ => None } }
    #[wasm_bindgen(getter)] pub fn branch(&self) -> Option<String> { match *self.inner { ck::CanonicalRankError::ProtocolDebt { branch, .. } => Some(branch.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn reason(&self) -> Option<String> { match *self.inner { ck::CanonicalRankError::ProtocolDebt { reason, .. } => Some(reason.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn carrier(&self) -> Option<usize> { match *self.inner { ck::CanonicalRankError::TemplateAttachmentRemap { carrier, .. } => Some(carrier.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn rings(&self) -> Option<CanonicalRingsPayload> { match self.inner.as_ref() { ck::CanonicalRankError::RingFinding(value) => Some(ring_errors::RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<CanonicalValencePayload> { match self.inner.as_ref() { ck::CanonicalRankError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = CanonicalCause, cause = cause);
