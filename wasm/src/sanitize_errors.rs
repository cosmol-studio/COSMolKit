//! Sanitization stage and exact contextual error causes.
pub(crate) use crate::enums::SanitizeStage;
use crate::{
    aromaticity_errors::AromaticityError,
    atropisomer_errors::AtropisomerError,
    errors::{enum_error, exception},
    kekulize_errors::KekulizeError,
    model_errors::{
        AtomPropertyError, BondValueError, MoleculePropertyError, TopologyValidationError,
    },
    query_errors::QueryStateError,
    ring_errors::RingFindingError,
    sanitize_leaf_errors::{AdjustHsError, HybridizationError, RadicalError},
    sanitize_support_errors::{CleanupError, ConjugationError, PropertyCacheError},
    stereo_errors::StereoError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "MoleculePropertyError | BondValueError | AtomPropertyError | TopologyValidationError | QueryStateError | CleanupError | PropertyCacheError | RingFindingError | KekulizeError | RadicalError | AromaticityError | ConjugationError | HybridizationError | AtropisomerError | StereoError | AdjustHsError"
    )]
    #[derive(Clone)]
    pub type SanitizeCause;
}
pub(crate) fn cause(value: &ck::SanitizeError) -> Option<SanitizeCause> {
    match value {
        ck::SanitizeError::MoleculeProperty(value) => {
            Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::SanitizeError::BondProperty(value) => {
            Some(BondValueError::from_rust(value.clone()).unchecked_into())
        }
        ck::SanitizeError::AtomProperty(value) => {
            Some(AtomPropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::SanitizeError::InvalidTopology { source, .. } => {
            Some(TopologyValidationError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::InvalidQueryState { source, .. } => {
            Some(QueryStateError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Cleanup { source, .. } => {
            Some(CleanupError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Properties { source, .. } => {
            Some(PropertyCacheError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Rings { source, .. } => {
            Some(RingFindingError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Kekulize { source, .. } => {
            Some(KekulizeError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Radicals { source, .. } => {
            Some(RadicalError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Aromaticity { source, .. } => {
            Some(AromaticityError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Conjugation { source, .. } => {
            Some(ConjugationError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Hybridization { source, .. } => {
            Some(HybridizationError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Atropisomers { source, .. } => {
            Some(AtropisomerError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::Chirality { source, .. } => {
            Some(StereoError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::AdjustHs { source, .. } => {
            Some(AdjustHsError::from_rust(source.clone()).unchecked_into())
        }
        ck::SanitizeError::InvalidOperations { .. } => None,
    }
}
enum_error!(SanitizeError, SanitizeErrorKind, ck::SanitizeError, {
    MoleculeProperty => ck::SanitizeError::MoleculeProperty(_), BondProperty => ck::SanitizeError::BondProperty(_), AtomProperty => ck::SanitizeError::AtomProperty(_),
    InvalidOperations => ck::SanitizeError::InvalidOperations { .. }, InvalidTopology => ck::SanitizeError::InvalidTopology { .. }, InvalidQueryState => ck::SanitizeError::InvalidQueryState { .. },
    Cleanup => ck::SanitizeError::Cleanup { .. }, Properties => ck::SanitizeError::Properties { .. }, Rings => ck::SanitizeError::Rings { .. }, Kekulize => ck::SanitizeError::Kekulize { .. },
    Radicals => ck::SanitizeError::Radicals { .. }, Aromaticity => ck::SanitizeError::Aromaticity { .. }, Conjugation => ck::SanitizeError::Conjugation { .. },
    Hybridization => ck::SanitizeError::Hybridization { .. }, Atropisomers => ck::SanitizeError::Atropisomers { .. }, Chirality => ck::SanitizeError::Chirality { .. }, AdjustHs => ck::SanitizeError::AdjustHs { .. },
}, {
    #[wasm_bindgen(getter)] pub fn bits(&self) -> Option<u32> { match *self.inner { ck::SanitizeError::InvalidOperations { bits, .. } => Some(bits), _ => None } }
    #[wasm_bindgen(getter, js_name = unknownBits)] pub fn unknown_bits(&self) -> Option<u32> { match *self.inner { ck::SanitizeError::InvalidOperations { unknown_bits, .. } => Some(unknown_bits), _ => None } }
    #[wasm_bindgen(getter)] pub fn stage(&self) -> Option<SanitizeStage> { match *self.inner {
        ck::SanitizeError::InvalidTopology { stage, .. } | ck::SanitizeError::InvalidQueryState { stage, .. } | ck::SanitizeError::Cleanup { stage, .. } | ck::SanitizeError::Properties { stage, .. }
        | ck::SanitizeError::Rings { stage, .. } | ck::SanitizeError::Kekulize { stage, .. } | ck::SanitizeError::Radicals { stage, .. } | ck::SanitizeError::Aromaticity { stage, .. }
        | ck::SanitizeError::Conjugation { stage, .. } | ck::SanitizeError::Hybridization { stage, .. } | ck::SanitizeError::Atropisomers { stage, .. } | ck::SanitizeError::Chirality { stage, .. } | ck::SanitizeError::AdjustHs { stage, .. } => Some(stage.into()), _ => None,
    } }
}, source = SanitizeCause, cause = cause);
