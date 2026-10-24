//! Kekulization error values and native transparent-source distinctions.
use crate::{
    canonical_errors,
    errors::{enum_error, exception},
    model_errors::TopologyValidationError,
    query_errors, ring_errors, valence_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | MappingValidationError | MoleculePropertyError | StereoGroupError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type KekulizeCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type KekulizeTopology;
    #[wasm_bindgen(typescript_type = "QueryStateError")]
    pub type KekulizeQuery;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type KekulizeRings;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type KekulizeValence;
    #[wasm_bindgen(typescript_type = "CanonicalRankError")]
    pub type KekulizeRank;
}
pub(crate) fn cause(value: &ck::KekulizeError) -> Option<KekulizeCause> {
    match value {
        ck::KekulizeError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::KekulizeError::InvalidQueryState(value) => {
            query_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::KekulizeError::RingFinding(value) => {
            ring_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::KekulizeError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::KekulizeError::CanonicalRank(value) => {
            canonical_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::KekulizeError::AtomSelectionLength { .. }
        | ck::KekulizeError::BondSelectionLength { .. }
        | ck::KekulizeError::CandidateAtomOutOfRange { .. }
        | ck::KekulizeError::DuplicateCandidateAtom { .. }
        | ck::KekulizeError::MatchingStateLength { .. }
        | ck::KekulizeError::DuplicateDoneAtom { .. }
        | ck::KekulizeError::MissingCandidateBond { .. }
        | ck::KekulizeError::MissingBacktrackAnchor { .. }
        | ck::KekulizeError::QuestionSubsetOverflow { .. }
        | ck::KekulizeError::AromaticAtomOutsideRing { .. }
        | ck::KekulizeError::NotKekulizable { .. }
        | ck::KekulizeError::PostconditionValenceMismatch { .. }
        | ck::KekulizeError::UnsupportedQueryState { .. }
        | ck::KekulizeError::IntegerOverflow { .. } => None,
    }
}
enum_error!(KekulizeError, KekulizeErrorKind, ck::KekulizeError, {
    InvalidTopology => ck::KekulizeError::InvalidTopology(_), AtomSelectionLength => ck::KekulizeError::AtomSelectionLength { .. },
    BondSelectionLength => ck::KekulizeError::BondSelectionLength { .. }, CandidateAtomOutOfRange => ck::KekulizeError::CandidateAtomOutOfRange { .. },
    DuplicateCandidateAtom => ck::KekulizeError::DuplicateCandidateAtom { .. }, MatchingStateLength => ck::KekulizeError::MatchingStateLength { .. },
    DuplicateDoneAtom => ck::KekulizeError::DuplicateDoneAtom { .. }, MissingCandidateBond => ck::KekulizeError::MissingCandidateBond { .. },
    MissingBacktrackAnchor => ck::KekulizeError::MissingBacktrackAnchor { .. }, QuestionSubsetOverflow => ck::KekulizeError::QuestionSubsetOverflow { .. },
    AromaticAtomOutsideRing => ck::KekulizeError::AromaticAtomOutsideRing { .. }, NotKekulizable => ck::KekulizeError::NotKekulizable { .. },
    PostconditionValenceMismatch => ck::KekulizeError::PostconditionValenceMismatch { .. }, UnsupportedQueryState => ck::KekulizeError::UnsupportedQueryState { .. },
    InvalidQueryState => ck::KekulizeError::InvalidQueryState(_), IntegerOverflow => ck::KekulizeError::IntegerOverflow { .. },
    RingFinding => ck::KekulizeError::RingFinding(_), Valence => ck::KekulizeError::Valence(_), CanonicalRank => ck::KekulizeError::CanonicalRank(_),
}, {
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::KekulizeError::AtomSelectionLength { expected, .. } | ck::KekulizeError::BondSelectionLength { expected, .. } | ck::KekulizeError::MatchingStateLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::KekulizeError::AtomSelectionLength { actual, .. } | ck::KekulizeError::BondSelectionLength { actual, .. } | ck::KekulizeError::MatchingStateLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::KekulizeError::CandidateAtomOutOfRange { atom, .. } | ck::KekulizeError::DuplicateCandidateAtom { atom } | ck::KekulizeError::DuplicateDoneAtom { atom }
        | ck::KekulizeError::MissingBacktrackAnchor { atom } | ck::KekulizeError::AromaticAtomOutsideRing { atom } | ck::KekulizeError::PostconditionValenceMismatch { atom, .. }
        | ck::KekulizeError::IntegerOverflow { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::KekulizeError::CandidateAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::KekulizeError::MatchingStateLength { field, .. } | ck::KekulizeError::IntegerOverflow { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn begin(&self) -> Option<usize> { match *self.inner { ck::KekulizeError::MissingCandidateBond { begin, .. } => Some(begin.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn end(&self) -> Option<usize> { match *self.inner { ck::KekulizeError::MissingCandidateBond { end, .. } => Some(end.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn questions(&self) -> Option<usize> { match *self.inner { ck::KekulizeError::QuestionSubsetOverflow { questions, .. } => Some(questions), _ => None } }
    #[wasm_bindgen(getter, js_name = bitWidth)] pub fn bit_width(&self) -> Option<u32> { match *self.inner { ck::KekulizeError::QuestionSubsetOverflow { bit_width, .. } => Some(bit_width), _ => None } }
    #[wasm_bindgen(getter, js_name = problemAtoms)] pub fn problem_atoms(&self) -> Option<Vec<u32>> { match self.inner.as_ref() { ck::KekulizeError::NotKekulizable { problem_atoms } => Some(problem_atoms.iter().map(|id| id.index() as u32).collect()), _ => None } }
    #[wasm_bindgen(getter)] pub fn before(&self) -> Option<i32> { match *self.inner { ck::KekulizeError::PostconditionValenceMismatch { before, .. } => Some(before), _ => None } }
    #[wasm_bindgen(getter)] pub fn after(&self) -> Option<i32> { match *self.inner { ck::KekulizeError::PostconditionValenceMismatch { after, .. } => Some(after), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::KekulizeError::UnsupportedQueryState { bond, .. } => Some(bond.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn detail(&self) -> Option<String> { match *self.inner { ck::KekulizeError::UnsupportedQueryState { detail, .. } => Some(detail.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<KekulizeTopology> { match self.inner.as_ref() { ck::KekulizeError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = queryState)] pub fn query_state(&self) -> Option<KekulizeQuery> { match self.inner.as_ref() { ck::KekulizeError::InvalidQueryState(value) => Some(query_errors::QueryStateError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = ringFinding)] pub fn ring_finding(&self) -> Option<KekulizeRings> { match self.inner.as_ref() { ck::KekulizeError::RingFinding(value) => Some(ring_errors::RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<KekulizeValence> { match self.inner.as_ref() { ck::KekulizeError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = canonicalRank)] pub fn canonical_rank(&self) -> Option<KekulizeRank> { match self.inner.as_ref() { ck::KekulizeError::CanonicalRank(value) => Some(canonical_errors::CanonicalRankError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = KekulizeCause, cause = cause);
