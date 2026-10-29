//! Tetrahedral order failures preserve native fields and actual sources.
use crate::{
    enums::ChiralTag,
    errors::{enum_error, exception},
    model_errors::{AtomPropertyError, MappingValidationError, TopologyValidationError},
    property_errors::PropertyUIntReadError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "PropertyUIntReadError | AtomPropertyError")]
    #[derive(Clone)]
    pub type StereoOrderCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type StereoOrderTopology;
    #[wasm_bindgen(typescript_type = "MappingValidationError")]
    pub type StereoOrderMapping;
}
pub(crate) fn cause(value: &ck::StereoOrderError) -> Option<StereoOrderCause> {
    match value {
        ck::StereoOrderError::ChiralPermutationRead(value) => {
            Some(PropertyUIntReadError::from_rust(value.clone()).unchecked_into())
        }
        ck::StereoOrderError::ChiralPermutationWrite(value) => {
            Some(AtomPropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::StereoOrderError::InvalidTopology(_)
        | ck::StereoOrderError::InvalidMapping(_)
        | ck::StereoOrderError::ConflictingChiralPermutation { .. }
        | ck::StereoOrderError::PermutationLength { .. }
        | ck::StereoOrderError::MissingProbeValue { .. }
        | ck::StereoOrderError::BondIndexSourceWidth { .. }
        | ck::StereoOrderError::CenterOutOfRange { .. }
        | ck::StereoOrderError::BondOutOfRange { .. }
        | ck::StereoOrderError::BondNotIncident { .. }
        | ck::StereoOrderError::TooManyLigands { .. }
        | ck::StereoOrderError::DuplicateBondLigand { .. }
        | ck::StereoOrderError::MultipleImplicitLigands { .. }
        | ck::StereoOrderError::UnsupportedTetrahedralTag { .. }
        | ck::StereoOrderError::RemovedCenter { .. }
        | ck::StereoOrderError::RemovedLigand { .. }
        | ck::StereoOrderError::ImplicitReplacementOutOfRange { .. }
        | ck::StereoOrderError::ImplicitReplacementRetained { .. }
        | ck::StereoOrderError::ImplicitReplacementNotInSource { .. } => None,
    }
}
enum_error!(StereoOrderError, StereoOrderErrorKind, ck::StereoOrderError, {
    ChiralPermutationRead => ck::StereoOrderError::ChiralPermutationRead(_), ChiralPermutationWrite => ck::StereoOrderError::ChiralPermutationWrite(_),
    ConflictingChiralPermutation => ck::StereoOrderError::ConflictingChiralPermutation { .. }, PermutationLength => ck::StereoOrderError::PermutationLength { .. },
    MissingProbeValue => ck::StereoOrderError::MissingProbeValue { .. }, BondIndexSourceWidth => ck::StereoOrderError::BondIndexSourceWidth { .. },
    InvalidTopology => ck::StereoOrderError::InvalidTopology(_), InvalidMapping => ck::StereoOrderError::InvalidMapping(_),
    CenterOutOfRange => ck::StereoOrderError::CenterOutOfRange { .. }, BondOutOfRange => ck::StereoOrderError::BondOutOfRange { .. },
    BondNotIncident => ck::StereoOrderError::BondNotIncident { .. }, TooManyLigands => ck::StereoOrderError::TooManyLigands { .. },
    DuplicateBondLigand => ck::StereoOrderError::DuplicateBondLigand { .. }, MultipleImplicitLigands => ck::StereoOrderError::MultipleImplicitLigands { .. },
    UnsupportedTetrahedralTag => ck::StereoOrderError::UnsupportedTetrahedralTag { .. }, RemovedCenter => ck::StereoOrderError::RemovedCenter { .. },
    RemovedLigand => ck::StereoOrderError::RemovedLigand { .. }, ImplicitReplacementOutOfRange => ck::StereoOrderError::ImplicitReplacementOutOfRange { .. },
    ImplicitReplacementRetained => ck::StereoOrderError::ImplicitReplacementRetained { .. }, ImplicitReplacementNotInSource => ck::StereoOrderError::ImplicitReplacementNotInSource { .. },
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::ConflictingChiralPermutation { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn typed(&self) -> Option<u32> { match *self.inner { ck::StereoOrderError::ConflictingChiralPermutation { typed, .. } => Some(typed), _ => None } }
    #[wasm_bindgen(getter)] pub fn raw(&self) -> Option<u32> { match *self.inner { ck::StereoOrderError::ConflictingChiralPermutation { raw, .. } => Some(raw), _ => None } }
    #[wasm_bindgen(getter)] pub fn reference(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::PermutationLength { reference, .. } => Some(reference), _ => None } }
    #[wasm_bindgen(getter)] pub fn probe(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::PermutationLength { probe, .. } => Some(probe), _ => None } }
    #[wasm_bindgen(getter, js_name = referencePosition)] pub fn reference_position(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::MissingProbeValue { reference_position } => Some(reference_position), _ => None } }
    #[wasm_bindgen(getter, js_name = bondIndex)] pub fn bond_index(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::BondIndexSourceWidth { bond_index } => Some(bond_index), _ => None } }
    #[wasm_bindgen(getter)] pub fn center(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::CenterOutOfRange { center, .. } | ck::StereoOrderError::BondNotIncident { center, .. } | ck::StereoOrderError::RemovedCenter { center } => Some(center.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::StereoOrderError::BondOutOfRange { bond, .. } | ck::StereoOrderError::BondNotIncident { bond, .. } | ck::StereoOrderError::DuplicateBondLigand { bond, .. }
        | ck::StereoOrderError::RemovedLigand { bond } | ck::StereoOrderError::ImplicitReplacementOutOfRange { bond, .. } | ck::StereoOrderError::ImplicitReplacementRetained { bond, .. }
        | ck::StereoOrderError::ImplicitReplacementNotInSource { bond } => Some(bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::CenterOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::BondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::TooManyLigands { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn maximum(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::TooManyLigands { maximum, .. } => Some(maximum), _ => None } }
    #[wasm_bindgen(getter)] pub fn position(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::DuplicateBondLigand { position, .. } => Some(position), _ => None } }
    #[wasm_bindgen(getter)] pub fn first(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::MultipleImplicitLigands { first, .. } => Some(first), _ => None } }
    #[wasm_bindgen(getter)] pub fn second(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::MultipleImplicitLigands { second, .. } => Some(second), _ => None } }
    #[wasm_bindgen(getter)] pub fn tag(&self) -> Option<ChiralTag> { match *self.inner { ck::StereoOrderError::UnsupportedTetrahedralTag { tag } => Some(tag.into()), _ => None } }
    #[wasm_bindgen(getter, js_name = oldBondCount)] pub fn old_bond_count(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::ImplicitReplacementOutOfRange { old_bond_count, .. } => Some(old_bond_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn mapped(&self) -> Option<usize> { match *self.inner { ck::StereoOrderError::ImplicitReplacementRetained { mapped, .. } => Some(mapped.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<StereoOrderTopology> { match self.inner.as_ref() { ck::StereoOrderError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn mapping(&self) -> Option<StereoOrderMapping> { match self.inner.as_ref() { ck::StereoOrderError::InvalidMapping(value) => Some(MappingValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = StereoOrderCause, cause = cause);
