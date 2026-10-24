//! Aromaticity diagnostics use the native enum, including transparent payloads.
use crate::enums::AromaticityModel;
use crate::{
    errors::{enum_error, exception},
    kekulize_errors,
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
    pub type AromaticityCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type AromaticityTopology;
    #[wasm_bindgen(typescript_type = "QueryStateError")]
    pub type AromaticityQuery;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type AromaticityRings;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type AromaticityValence;
    #[wasm_bindgen(typescript_type = "KekulizeError")]
    pub type AromaticityKekulize;
}
pub(crate) fn cause(value: &ck::AromaticityError) -> Option<AromaticityCause> {
    match value {
        ck::AromaticityError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::AromaticityError::InvalidQueryState(value) => {
            query_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::AromaticityError::RingFinding(value) => {
            ring_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::AromaticityError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::AromaticityError::Kekulize(value) => {
            kekulize_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::AromaticityError::RingInfoNotInitialized
        | ck::AromaticityError::RingInfoDimensionMismatch { .. }
        | ck::AromaticityError::RingTableLengthMismatch { .. }
        | ck::AromaticityError::RingRowLengthMismatch { .. }
        | ck::AromaticityError::RingAtomOutOfRange { .. }
        | ck::AromaticityError::RingBondOutOfRange { .. }
        | ck::AromaticityError::AtomOutOfRange { .. }
        | ck::AromaticityError::BondOutOfRange { .. }
        | ck::AromaticityError::ValenceAssignmentLength { .. }
        | ck::AromaticityError::InvalidValenceRow { .. }
        | ck::AromaticityError::ExpectedRingBondNotFound { .. }
        | ck::AromaticityError::UnsupportedModel { .. }
        | ck::AromaticityError::UnsupportedState { .. }
        | ck::AromaticityError::IntegerOverflow { .. }
        | ck::AromaticityError::UnexpectedTopologyShape { .. }
        | ck::AromaticityError::AtomIdentityChanged { .. }
        | ck::AromaticityError::BondIdentityChanged { .. }
        | ck::AromaticityError::AromaticRingCountOutOfRange { .. } => None,
    }
}
enum_error!(AromaticityError, AromaticityErrorKind, ck::AromaticityError, {
    InvalidTopology => ck::AromaticityError::InvalidTopology(_), InvalidQueryState => ck::AromaticityError::InvalidQueryState(_), RingInfoNotInitialized => ck::AromaticityError::RingInfoNotInitialized,
    RingInfoDimensionMismatch => ck::AromaticityError::RingInfoDimensionMismatch { .. }, RingTableLengthMismatch => ck::AromaticityError::RingTableLengthMismatch { .. },
    RingRowLengthMismatch => ck::AromaticityError::RingRowLengthMismatch { .. }, RingAtomOutOfRange => ck::AromaticityError::RingAtomOutOfRange { .. }, RingBondOutOfRange => ck::AromaticityError::RingBondOutOfRange { .. },
    AtomOutOfRange => ck::AromaticityError::AtomOutOfRange { .. }, BondOutOfRange => ck::AromaticityError::BondOutOfRange { .. }, ValenceAssignmentLength => ck::AromaticityError::ValenceAssignmentLength { .. },
    InvalidValenceRow => ck::AromaticityError::InvalidValenceRow { .. }, ExpectedRingBondNotFound => ck::AromaticityError::ExpectedRingBondNotFound { .. }, UnsupportedModel => ck::AromaticityError::UnsupportedModel { .. },
    UnsupportedState => ck::AromaticityError::UnsupportedState { .. }, IntegerOverflow => ck::AromaticityError::IntegerOverflow { .. }, UnexpectedTopologyShape => ck::AromaticityError::UnexpectedTopologyShape { .. },
    AtomIdentityChanged => ck::AromaticityError::AtomIdentityChanged { .. }, BondIdentityChanged => ck::AromaticityError::BondIdentityChanged { .. }, AromaticRingCountOutOfRange => ck::AromaticityError::AromaticRingCountOutOfRange { .. },
    RingFinding => ck::AromaticityError::RingFinding(_), Valence => ck::AromaticityError::Valence(_), Kekulize => ck::AromaticityError::Kekulize(_),
}, {
    #[wasm_bindgen(getter, js_name = ringAtomCount)] pub fn ring_atom_count(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingInfoDimensionMismatch { ring_atom_count, .. } => Some(ring_atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = ringBondCount)] pub fn ring_bond_count(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingInfoDimensionMismatch { ring_bond_count, .. } => Some(ring_bond_count), _ => None } }
    #[wasm_bindgen(getter, js_name = topologyAtomCount)] pub fn topology_atom_count(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingInfoDimensionMismatch { topology_atom_count, .. } => Some(topology_atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = topologyBondCount)] pub fn topology_bond_count(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingInfoDimensionMismatch { topology_bond_count, .. } => Some(topology_bond_count), _ => None } }
    #[wasm_bindgen(getter, js_name = atomRows)] pub fn atom_rows(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingTableLengthMismatch { atom_rows, .. } => Some(atom_rows), _ => None } }
    #[wasm_bindgen(getter, js_name = bondRows)] pub fn bond_rows(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingTableLengthMismatch { bond_rows, .. } => Some(bond_rows), _ => None } }
    #[wasm_bindgen(getter, js_name = ringIndex)] pub fn ring_index(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingRowLengthMismatch { ring_index, .. } | ck::AromaticityError::RingAtomOutOfRange { ring_index, .. } | ck::AromaticityError::RingBondOutOfRange { ring_index, .. } => Some(ring_index), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingRowLengthMismatch { atom_count, .. } | ck::AromaticityError::RingAtomOutOfRange { atom_count, .. } | ck::AromaticityError::AtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingRowLengthMismatch { bond_count, .. } | ck::AromaticityError::RingBondOutOfRange { bond_count, .. } | ck::AromaticityError::BondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingAtomOutOfRange { atom, .. } | ck::AromaticityError::AtomOutOfRange { atom, .. } | ck::AromaticityError::InvalidValenceRow { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::RingBondOutOfRange { bond, .. } | ck::AromaticityError::BondOutOfRange { bond, .. } => Some(bond.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::AromaticityError::ValenceAssignmentLength { field, .. } | ck::AromaticityError::InvalidValenceRow { field, .. } | ck::AromaticityError::UnsupportedState { field, .. } | ck::AromaticityError::IntegerOverflow { field } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::ValenceAssignmentLength { actual, .. } | ck::AromaticityError::AromaticRingCountOutOfRange { actual, .. } => Some(actual), ck::AromaticityError::AtomIdentityChanged { actual, .. } => Some(actual.index()), ck::AromaticityError::BondIdentityChanged { actual, .. } => Some(actual.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::ValenceAssignmentLength { expected, .. } => Some(expected), ck::AromaticityError::AtomIdentityChanged { expected, .. } => Some(expected.index()), ck::AromaticityError::BondIdentityChanged { expected, .. } => Some(expected.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<i32> { match *self.inner { ck::AromaticityError::InvalidValenceRow { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn begin(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::ExpectedRingBondNotFound { begin, .. } => Some(begin.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn end(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::ExpectedRingBondNotFound { end, .. } => Some(end.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn model(&self) -> Option<AromaticityModel> { match *self.inner { ck::AromaticityError::UnsupportedModel { model, .. } | ck::AromaticityError::AromaticRingCountOutOfRange { model, .. } => Some(model.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn detail(&self) -> Option<String> { match *self.inner { ck::AromaticityError::UnsupportedModel { detail, .. } | ck::AromaticityError::UnsupportedState { detail, .. } => Some(detail.into()), _ => None } }
    #[wasm_bindgen(getter, js_name = inputAtoms)] pub fn input_atoms(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::UnexpectedTopologyShape { input_atoms, .. } => Some(input_atoms), _ => None } }
    #[wasm_bindgen(getter, js_name = inputBonds)] pub fn input_bonds(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::UnexpectedTopologyShape { input_bonds, .. } => Some(input_bonds), _ => None } }
    #[wasm_bindgen(getter, js_name = outputAtoms)] pub fn output_atoms(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::UnexpectedTopologyShape { output_atoms, .. } => Some(output_atoms), _ => None } }
    #[wasm_bindgen(getter, js_name = outputBonds)] pub fn output_bonds(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::UnexpectedTopologyShape { output_bonds, .. } => Some(output_bonds), _ => None } }
    #[wasm_bindgen(getter)] pub fn row(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::AtomIdentityChanged { row, .. } | ck::AromaticityError::BondIdentityChanged { row, .. } => Some(row), _ => None } }
    #[wasm_bindgen(getter)] pub fn maximum(&self) -> Option<usize> { match *self.inner { ck::AromaticityError::AromaticRingCountOutOfRange { maximum, .. } => Some(maximum), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<AromaticityTopology> { match self.inner.as_ref() { ck::AromaticityError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = queryState)] pub fn query_state(&self) -> Option<AromaticityQuery> { match self.inner.as_ref() { ck::AromaticityError::InvalidQueryState(value) => Some(query_errors::QueryStateError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = ringFinding)] pub fn ring_finding(&self) -> Option<AromaticityRings> { match self.inner.as_ref() { ck::AromaticityError::RingFinding(value) => Some(ring_errors::RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<AromaticityValence> { match self.inner.as_ref() { ck::AromaticityError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn kekulize(&self) -> Option<AromaticityKekulize> { match self.inner.as_ref() { ck::AromaticityError::Kekulize(value) => Some(kekulize_errors::KekulizeError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = AromaticityCause, cause = cause);
