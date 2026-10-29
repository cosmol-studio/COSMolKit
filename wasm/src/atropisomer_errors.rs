//! Axial stereochemistry errors preserve native payload and cause semantics.
use crate::{
    enums::BondStereo,
    errors::{enum_error, exception},
    model_errors::{
        BondValueError, CoordinateValidationError, StereoGroupError, TopologyValidationError,
    },
    ring_errors::RingFindingError,
    sanitize_leaf_errors::HybridizationError,
    sanitize_support_errors::ConjugationError,
    stereo_errors::StereoError,
    valence_errors::ValenceError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | ValenceError | ConjugationError | HybridizationError | RingFindingError | TopologyValidationError | CoordinateValidationError | BondValueError"
    )]
    #[derive(Clone)]
    pub type AtropisomerCause;
    #[wasm_bindgen(typescript_type = "StereoError")]
    pub type AtropisomerNormalization;
}
pub(crate) fn cause(value: &ck::AtropisomerError) -> Option<AtropisomerCause> {
    match value {
        ck::AtropisomerError::StereoGroup(value) => {
            Some(StereoGroupError::from_rust(*value).unchecked_into())
        }
        ck::AtropisomerError::Valence(value) => {
            Some(ValenceError::from_rust(value.clone()).unchecked_into())
        }
        ck::AtropisomerError::Conjugation(value) => {
            Some(ConjugationError::from_rust(value.clone()).unchecked_into())
        }
        ck::AtropisomerError::Hybridization(value) => {
            Some(HybridizationError::from_rust(value.clone()).unchecked_into())
        }
        ck::AtropisomerError::RingFinding(value) => {
            Some(RingFindingError::from_rust(value.clone()).unchecked_into())
        }
        ck::AtropisomerError::InvalidTopology { source } => {
            Some(TopologyValidationError::from_rust(source.clone()).unchecked_into())
        }
        ck::AtropisomerError::InvalidCoordinates { source } => {
            Some(CoordinateValidationError::from_rust(source.clone()).unchecked_into())
        }
        ck::AtropisomerError::BondUpdate { source, .. } => {
            Some(BondValueError::from_rust(source.clone()).unchecked_into())
        }
        ck::AtropisomerError::Normalization(_)
        | ck::AtropisomerError::SourcePrecondition { .. }
        | ck::AtropisomerError::CandidateTraversalOutOfRange { .. }
        | ck::AtropisomerError::CandidateTraversalDuplicate { .. }
        | ck::AtropisomerError::CandidateTraversalUnexpected { .. }
        | ck::AtropisomerError::CandidateTraversalMissing { .. }
        | ck::AtropisomerError::CarrierCount { .. }
        | ck::AtropisomerError::ConformerNotThreeDimensional { .. }
        | ck::AtropisomerError::RingInfoNotSssr
        | ck::AtropisomerError::RingAtomOutOfRange { .. }
        | ck::AtropisomerError::RingBondOutOfRange { .. }
        | ck::AtropisomerError::RingAtomRowCount { .. }
        | ck::AtropisomerError::RingBondRowCount { .. }
        | ck::AtropisomerError::HybridizationAssignmentLength { .. }
        | ck::AtropisomerError::StereoGroupAtomOutOfRange { .. }
        | ck::AtropisomerError::StereoGroupBondOutOfRange { .. }
        | ck::AtropisomerError::AssignmentBondOutOfRange { .. }
        | ck::AtropisomerError::AxialBondOutOfRange { .. }
        | ck::AtropisomerError::InvalidAtropisomerStereo { .. } => None,
    }
}
enum_error!(AtropisomerError, AtropisomerErrorKind, ck::AtropisomerError, {
    StereoGroup => ck::AtropisomerError::StereoGroup(_), SourcePrecondition => ck::AtropisomerError::SourcePrecondition { .. },
    CandidateTraversalOutOfRange => ck::AtropisomerError::CandidateTraversalOutOfRange { .. }, CandidateTraversalDuplicate => ck::AtropisomerError::CandidateTraversalDuplicate { .. },
    CandidateTraversalUnexpected => ck::AtropisomerError::CandidateTraversalUnexpected { .. }, CandidateTraversalMissing => ck::AtropisomerError::CandidateTraversalMissing { .. },
    Valence => ck::AtropisomerError::Valence(_), Conjugation => ck::AtropisomerError::Conjugation(_), Hybridization => ck::AtropisomerError::Hybridization(_),
    Normalization => ck::AtropisomerError::Normalization(_), CarrierCount => ck::AtropisomerError::CarrierCount { .. }, RingFinding => ck::AtropisomerError::RingFinding(_),
    InvalidTopology => ck::AtropisomerError::InvalidTopology { .. }, InvalidCoordinates => ck::AtropisomerError::InvalidCoordinates { .. },
    ConformerNotThreeDimensional => ck::AtropisomerError::ConformerNotThreeDimensional { .. }, RingInfoNotSssr => ck::AtropisomerError::RingInfoNotSssr,
    RingAtomOutOfRange => ck::AtropisomerError::RingAtomOutOfRange { .. }, RingBondOutOfRange => ck::AtropisomerError::RingBondOutOfRange { .. },
    RingAtomRowCount => ck::AtropisomerError::RingAtomRowCount { .. }, RingBondRowCount => ck::AtropisomerError::RingBondRowCount { .. },
    HybridizationAssignmentLength => ck::AtropisomerError::HybridizationAssignmentLength { .. }, BondUpdate => ck::AtropisomerError::BondUpdate { .. },
    StereoGroupAtomOutOfRange => ck::AtropisomerError::StereoGroupAtomOutOfRange { .. }, StereoGroupBondOutOfRange => ck::AtropisomerError::StereoGroupBondOutOfRange { .. },
    AssignmentBondOutOfRange => ck::AtropisomerError::AssignmentBondOutOfRange { .. }, AxialBondOutOfRange => ck::AtropisomerError::AxialBondOutOfRange { .. },
    InvalidAtropisomerStereo => ck::AtropisomerError::InvalidAtropisomerStereo { .. },
}, {
    #[wasm_bindgen(getter)] pub fn message(&self) -> Option<String> { match *self.inner { ck::AtropisomerError::SourcePrecondition { message } => Some(message.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::AtropisomerError::CandidateTraversalOutOfRange { bond, .. } | ck::AtropisomerError::CandidateTraversalDuplicate { bond }
        | ck::AtropisomerError::CandidateTraversalUnexpected { bond } | ck::AtropisomerError::CandidateTraversalMissing { bond } | ck::AtropisomerError::RingBondOutOfRange { bond, .. }
        | ck::AtropisomerError::BondUpdate { bond, .. } | ck::AtropisomerError::StereoGroupBondOutOfRange { bond, .. } | ck::AtropisomerError::AssignmentBondOutOfRange { bond, .. }
        | ck::AtropisomerError::AxialBondOutOfRange { bond, .. } | ck::AtropisomerError::InvalidAtropisomerStereo { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::AtropisomerError::CarrierCount { atom, .. } | ck::AtropisomerError::RingAtomOutOfRange { atom, .. } | ck::AtropisomerError::StereoGroupAtomOutOfRange { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::AtropisomerError::CandidateTraversalOutOfRange { bond_count, .. } | ck::AtropisomerError::RingBondOutOfRange { bond_count, .. } | ck::AtropisomerError::StereoGroupBondOutOfRange { bond_count, .. } | ck::AtropisomerError::AssignmentBondOutOfRange { bond_count, .. } | ck::AtropisomerError::AxialBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::AtropisomerError::RingAtomOutOfRange { atom_count, .. } | ck::AtropisomerError::StereoGroupAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::AtropisomerError::RingAtomRowCount { actual, .. } | ck::AtropisomerError::RingBondRowCount { actual, .. } | ck::AtropisomerError::HybridizationAssignmentLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::AtropisomerError::RingAtomRowCount { expected, .. } | ck::AtropisomerError::RingBondRowCount { expected, .. } | ck::AtropisomerError::HybridizationAssignmentLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn count(&self) -> Option<usize> { match *self.inner { ck::AtropisomerError::CarrierCount { count, .. } => Some(count), _ => None } }
    #[wasm_bindgen(getter)] pub fn conformer(&self) -> Option<usize> { match *self.inner { ck::AtropisomerError::ConformerNotThreeDimensional { conformer } => Some(conformer), _ => None } }
    #[wasm_bindgen(getter)] pub fn stereo(&self) -> Option<BondStereo> { match *self.inner { ck::AtropisomerError::InvalidAtropisomerStereo { stereo, .. } => Some(stereo.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn normalization(&self) -> Option<AtropisomerNormalization> { match self.inner.as_ref() { ck::AtropisomerError::Normalization(value) => Some(StereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = AtropisomerCause, cause = cause);
