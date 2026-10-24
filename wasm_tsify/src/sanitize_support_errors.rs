//! Native cleanup/conjugation/property-cache error wrappers.
use crate::{
    aromaticity_errors, canonical_errors,
    errors::{enum_error, exception},
    model_errors::TopologyValidationError,
    valence_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | MappingValidationError | MoleculePropertyError | StereoGroupError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type ConjugationCause;
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | MoleculePropertyError | StereoGroupError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type CleanupCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type SupportTopology;
    #[wasm_bindgen(typescript_type = "CanonicalRankError")]
    pub type SupportRank;
    #[wasm_bindgen(typescript_type = "AromaticityError")]
    pub type SupportAromaticity;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    #[derive(Clone)]
    pub type SupportValence;
}
pub(crate) fn conjugation_cause(value: &ck::ConjugationError) -> Option<ConjugationCause> {
    match value {
        ck::ConjugationError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::ConjugationError::Aromaticity(value) => {
            aromaticity_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::ConjugationError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::ConjugationError::ValenceAssignmentLength { .. }
        | ck::ConjugationError::InvalidValenceRow { .. }
        | ck::ConjugationError::AtomOutOfRange { .. }
        | ck::ConjugationError::BondOutOfRange { .. }
        | ck::ConjugationError::IntegerOverflow { .. } => None,
    }
}
enum_error!(ConjugationError, ConjugationErrorKind, ck::ConjugationError, {
    InvalidTopology => ck::ConjugationError::InvalidTopology(_), ValenceAssignmentLength => ck::ConjugationError::ValenceAssignmentLength { .. },
    InvalidValenceRow => ck::ConjugationError::InvalidValenceRow { .. }, AtomOutOfRange => ck::ConjugationError::AtomOutOfRange { .. },
    BondOutOfRange => ck::ConjugationError::BondOutOfRange { .. }, IntegerOverflow => ck::ConjugationError::IntegerOverflow { .. },
    Aromaticity => ck::ConjugationError::Aromaticity(_), Valence => ck::ConjugationError::Valence(_),
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::ConjugationError::InvalidValenceRow { atom, .. } | ck::ConjugationError::AtomOutOfRange { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::ConjugationError::BondOutOfRange { bond, .. } => Some(bond.index()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::ConjugationError::AtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::ConjugationError::BondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::ConjugationError::ValenceAssignmentLength { field, .. } | ck::ConjugationError::InvalidValenceRow { field, .. } | ck::ConjugationError::IntegerOverflow { field } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::ConjugationError::ValenceAssignmentLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::ConjugationError::ValenceAssignmentLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<i32> { match *self.inner { ck::ConjugationError::InvalidValenceRow { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<SupportTopology> { match self.inner.as_ref() { ck::ConjugationError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn aromaticity(&self) -> Option<SupportAromaticity> { match self.inner.as_ref() { ck::ConjugationError::Aromaticity(value) => Some(aromaticity_errors::AromaticityError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<SupportValence> { match self.inner.as_ref() { ck::ConjugationError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = ConjugationCause, cause = conjugation_cause);

pub(crate) fn cleanup_cause(value: &ck::CleanupError) -> Option<CleanupCause> {
    match value {
        ck::CleanupError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::CleanupError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CleanupError::CanonicalRank(value) => {
            canonical_errors::cause(value).map(JsCast::unchecked_into)
        }
    }
}
enum_error!(CleanupError, CleanupErrorKind, ck::CleanupError, {
    InvalidTopology => ck::CleanupError::InvalidTopology(_), Valence => ck::CleanupError::Valence(_), CanonicalRank => ck::CleanupError::CanonicalRank(_),
}, {
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<SupportTopology> { match self.inner.as_ref() { ck::CleanupError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<SupportValence> { match self.inner.as_ref() { ck::CleanupError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = canonicalRank)] pub fn canonical_rank(&self) -> Option<SupportRank> { match self.inner.as_ref() { ck::CleanupError::CanonicalRank(value) => Some(canonical_errors::CanonicalRankError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = CleanupCause, cause = cleanup_cause);

enum_error!(PropertyCacheError, PropertyCacheErrorKind, ck::PropertyCacheError, {
    Valence => ck::PropertyCacheError::Valence(_),
}, {
    #[wasm_bindgen(getter)] pub fn valence(&self) -> SupportValence { match self.inner.as_ref() { ck::PropertyCacheError::Valence(value) => valence_errors::ValenceError::from_rust(value.clone()).unchecked_into() } }
}, source = SupportValence, cause = |value: &ck::PropertyCacheError| match value {
    ck::PropertyCacheError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()),
});
