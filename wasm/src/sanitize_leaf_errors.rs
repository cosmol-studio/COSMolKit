//! Sanitization leaf diagnostics share native topology-source projection.
use crate::{
    errors::{enum_error, exception},
    model_errors::TopologyValidationError,
    valence_errors::{self, ValenceError},
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    #[derive(Clone)]
    pub type SanitizeLeafCause;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type SanitizeLeafValence;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type HybridizationTopology;
}
pub(crate) fn adjust_cause(value: &ck::AdjustHsError) -> Option<SanitizeLeafCause> {
    match value {
        ck::AdjustHsError::InvalidTopology { source } => {
            Some(TopologyValidationError::from_rust(source.clone()).unchecked_into())
        }
        ck::AdjustHsError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::AdjustHsError::ValenceAssignmentLength { .. }
        | ck::AdjustHsError::InvalidValenceRow { .. }
        | ck::AdjustHsError::ExplicitHydrogenOverflow { .. } => None,
    }
}
enum_error!(AdjustHsError, AdjustHsErrorKind, ck::AdjustHsError, {
    InvalidTopology => ck::AdjustHsError::InvalidTopology { .. }, ValenceAssignmentLength => ck::AdjustHsError::ValenceAssignmentLength { .. },
    InvalidValenceRow => ck::AdjustHsError::InvalidValenceRow { .. }, ExplicitHydrogenOverflow => ck::AdjustHsError::ExplicitHydrogenOverflow { .. }, Valence => ck::AdjustHsError::Valence(_),
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::AdjustHsError::InvalidValenceRow { atom, .. } | ck::AdjustHsError::ExplicitHydrogenOverflow { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::AdjustHsError::ValenceAssignmentLength { field, .. } | ck::AdjustHsError::InvalidValenceRow { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::AdjustHsError::ValenceAssignmentLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::AdjustHsError::ValenceAssignmentLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<i32> { match *self.inner { ck::AdjustHsError::InvalidValenceRow { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter, js_name = originalExplicitHydrogens)] pub fn original_explicit_hydrogens(&self) -> Option<u8> { match *self.inner { ck::AdjustHsError::ExplicitHydrogenOverflow { original_explicit_hydrogens, .. } => Some(original_explicit_hydrogens), _ => None } }
    #[wasm_bindgen(getter, js_name = originalImplicitValence)] pub fn original_implicit_valence(&self) -> Option<i32> { match *self.inner { ck::AdjustHsError::ExplicitHydrogenOverflow { original_implicit_valence, .. } => Some(original_implicit_valence), _ => None } }
    #[wasm_bindgen(getter, js_name = recalculatedImplicitValence)] pub fn recalculated_implicit_valence(&self) -> Option<i32> { match *self.inner { ck::AdjustHsError::ExplicitHydrogenOverflow { recalculated_implicit_valence, .. } => Some(recalculated_implicit_valence), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<SanitizeLeafValence> { match self.inner.as_ref() { ck::AdjustHsError::Valence(value) => Some(ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = SanitizeLeafCause, cause = adjust_cause);

pub(crate) fn radical_cause(value: &ck::RadicalError) -> Option<SanitizeLeafCause> {
    match value {
        ck::RadicalError::InvalidTopology { source } => {
            Some(TopologyValidationError::from_rust(source.clone()).unchecked_into())
        }
        ck::RadicalError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::RadicalError::PeriodicTableLookup { .. }
        | ck::RadicalError::RadicalCountOutOfRange { .. } => None,
    }
}
enum_error!(RadicalError, RadicalErrorKind, ck::RadicalError, {
    InvalidTopology => ck::RadicalError::InvalidTopology { .. }, PeriodicTableLookup => ck::RadicalError::PeriodicTableLookup { .. }, Valence => ck::RadicalError::Valence(_), RadicalCountOutOfRange => ck::RadicalError::RadicalCountOutOfRange { .. },
}, {
    #[wasm_bindgen(getter, js_name = atomicNumber)] pub fn atomic_number(&self) -> Option<u8> { match *self.inner { ck::RadicalError::PeriodicTableLookup { atomic_number, .. } => Some(atomic_number), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::RadicalError::PeriodicTableLookup { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::RadicalError::RadicalCountOutOfRange { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn count(&self) -> Option<i32> { match *self.inner { ck::RadicalError::RadicalCountOutOfRange { count, .. } => Some(count), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<SanitizeLeafValence> { match self.inner.as_ref() { ck::RadicalError::Valence(value) => Some(ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = SanitizeLeafCause, cause = radical_cause);

pub(crate) fn hybridization_cause(value: &ck::HybridizationError) -> Option<SanitizeLeafCause> {
    match value {
        ck::HybridizationError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::HybridizationError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::HybridizationError::ConjugatedBondAssignmentLength { .. }
        | ck::HybridizationError::ValenceAssignmentLength { .. }
        | ck::HybridizationError::InvalidValenceRow { .. }
        | ck::HybridizationError::IntegerOverflow { .. } => None,
    }
}
enum_error!(HybridizationError, HybridizationErrorKind, ck::HybridizationError, {
    ConjugatedBondAssignmentLength => ck::HybridizationError::ConjugatedBondAssignmentLength { .. }, InvalidTopology => ck::HybridizationError::InvalidTopology(_),
    ValenceAssignmentLength => ck::HybridizationError::ValenceAssignmentLength { .. }, InvalidValenceRow => ck::HybridizationError::InvalidValenceRow { .. },
    IntegerOverflow => ck::HybridizationError::IntegerOverflow { .. }, Valence => ck::HybridizationError::Valence(_),
}, {
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::HybridizationError::ConjugatedBondAssignmentLength { actual, .. } | ck::HybridizationError::ValenceAssignmentLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::HybridizationError::ConjugatedBondAssignmentLength { expected, .. } | ck::HybridizationError::ValenceAssignmentLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::HybridizationError::ValenceAssignmentLength { field, .. } | ck::HybridizationError::InvalidValenceRow { field, .. } | ck::HybridizationError::IntegerOverflow { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::HybridizationError::InvalidValenceRow { atom, .. } | ck::HybridizationError::IntegerOverflow { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<i32> { match *self.inner { ck::HybridizationError::InvalidValenceRow { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<HybridizationTopology> { match self.inner.as_ref() { ck::HybridizationError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<SanitizeLeafValence> { match self.inner.as_ref() { ck::HybridizationError::Valence(value) => Some(ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = SanitizeLeafCause, cause = hybridization_cause);
