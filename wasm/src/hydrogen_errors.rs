//! Hydrogen errors retain contextual fields and source-free error payloads.
use crate::{
    enums::SdfPropertyListTarget,
    errors::{enum_error, exception},
    model_errors::{
        AtomPropertyError, BondValueError, CoordinateValidationError, MappingValidationError,
        MoleculePropertyError, TopologyValidationError,
    },
    property_errors,
    query_errors::QueryStateError,
    sanitize_errors::SanitizeError,
    stereo_order_errors::StereoOrderError,
    topology_errors::TopologyEditError,
    valence_errors::ValenceError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoOrderError | AtomPropertyError | BondValueError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type HydrogenCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type HydrogenTopology;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type HydrogenCoordinates;
    #[wasm_bindgen(typescript_type = "TopologyEditError")]
    pub type HydrogenEdit;
    #[wasm_bindgen(typescript_type = "MappingValidationError")]
    pub type HydrogenMapping;
    #[wasm_bindgen(typescript_type = "QueryStateError")]
    pub type HydrogenQuery;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type HydrogenValence;
    #[wasm_bindgen(typescript_type = "SanitizeError")]
    pub type HydrogenSanitize;
}
pub(crate) fn cause(value: &ck::HydrogenError) -> Option<HydrogenCause> {
    match value {
        ck::HydrogenError::StereoOrder(value) => {
            Some(StereoOrderError::from_rust(value.clone()).unchecked_into())
        }
        ck::HydrogenError::BondProperty(value) => {
            Some(BondValueError::from_rust(value.clone()).unchecked_into())
        }
        ck::HydrogenError::AtomProperty(value) => {
            Some(AtomPropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::HydrogenError::InvalidProperty(value) => {
            Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::HydrogenError::InvalidTopology(_)
        | ck::HydrogenError::InvalidCoordinates(_)
        | ck::HydrogenError::TopologyEdit(_)
        | ck::HydrogenError::InvalidMapping(_)
        | ck::HydrogenError::InvalidQueryState(_)
        | ck::HydrogenError::InvalidPropertyList { .. }
        | ck::HydrogenError::Unsupported { .. }
        | ck::HydrogenError::OnlyOnAtomOutOfRange { .. }
        | ck::HydrogenError::InvalidAdditionPlan { .. }
        | ck::HydrogenError::InvalidRemovalCandidate { .. }
        | ck::HydrogenError::MissingSourceState { .. }
        | ck::HydrogenError::ValenceAssignmentLength { .. }
        | ck::HydrogenError::ExplicitHydrogenOverflow { .. }
        | ck::HydrogenError::InvalidStereoTransition { .. }
        | ck::HydrogenError::CoordinatePlacement { .. }
        | ck::HydrogenError::Valence(_)
        | ck::HydrogenError::Sanitize(_) => None,
    }
}
enum_error!(HydrogenError, HydrogenErrorKind, ck::HydrogenError, {
    StereoOrder => ck::HydrogenError::StereoOrder(_), BondProperty => ck::HydrogenError::BondProperty(_), AtomProperty => ck::HydrogenError::AtomProperty(_),
    InvalidTopology => ck::HydrogenError::InvalidTopology(_), InvalidCoordinates => ck::HydrogenError::InvalidCoordinates(_), TopologyEdit => ck::HydrogenError::TopologyEdit(_),
    InvalidMapping => ck::HydrogenError::InvalidMapping(_), InvalidQueryState => ck::HydrogenError::InvalidQueryState(_), InvalidProperty => ck::HydrogenError::InvalidProperty(_),
    InvalidPropertyList => ck::HydrogenError::InvalidPropertyList { .. }, Unsupported => ck::HydrogenError::Unsupported { .. }, OnlyOnAtomOutOfRange => ck::HydrogenError::OnlyOnAtomOutOfRange { .. },
    InvalidAdditionPlan => ck::HydrogenError::InvalidAdditionPlan { .. }, InvalidRemovalCandidate => ck::HydrogenError::InvalidRemovalCandidate { .. },
    MissingSourceState => ck::HydrogenError::MissingSourceState { .. }, ValenceAssignmentLength => ck::HydrogenError::ValenceAssignmentLength { .. },
    ExplicitHydrogenOverflow => ck::HydrogenError::ExplicitHydrogenOverflow { .. }, InvalidStereoTransition => ck::HydrogenError::InvalidStereoTransition { .. },
    CoordinatePlacement => ck::HydrogenError::CoordinatePlacement { .. }, Valence => ck::HydrogenError::Valence(_), Sanitize => ck::HydrogenError::Sanitize(_),
}, {
    #[wasm_bindgen(getter)] pub fn name(&self) -> Result<Option<String>, JsValue> { match self.inner.as_ref() { ck::HydrogenError::InvalidPropertyList { name, .. } => property_errors::text(name).map(Some), _ => Ok(None) } }
    #[wasm_bindgen(getter)] pub fn target(&self) -> Option<SdfPropertyListTarget> { match *self.inner { ck::HydrogenError::InvalidPropertyList { target, .. } => Some(target.into()), _ => None } }
    #[wasm_bindgen(getter, js_name = expectedRows)] pub fn expected_rows(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::InvalidPropertyList { expected_rows, .. } => Some(expected_rows), _ => None } }
    #[wasm_bindgen(getter, js_name = actualRows)] pub fn actual_rows(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::InvalidPropertyList { actual_rows, .. } => Some(actual_rows), _ => None } }
    #[wasm_bindgen(getter)] pub fn operation(&self) -> Option<String> { match *self.inner { ck::HydrogenError::Unsupported { operation, .. } => Some(operation.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn reason(&self) -> Option<String> { match *self.inner {
        ck::HydrogenError::Unsupported { reason, .. } | ck::HydrogenError::InvalidAdditionPlan { reason, .. } | ck::HydrogenError::InvalidRemovalCandidate { reason, .. }
        | ck::HydrogenError::InvalidStereoTransition { reason, .. } | ck::HydrogenError::CoordinatePlacement { reason, .. } => Some(reason.into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::OnlyOnAtomOutOfRange { atom, .. } | ck::HydrogenError::InvalidRemovalCandidate { atom, .. } | ck::HydrogenError::ExplicitHydrogenOverflow { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::OnlyOnAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn addition(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::InvalidAdditionPlan { addition, .. } => addition, ck::HydrogenError::CoordinatePlacement { addition, .. } => Some(addition), _ => None } }
    #[wasm_bindgen(getter)] pub fn position(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::InvalidRemovalCandidate { position, .. } => Some(position), _ => None } }
    #[wasm_bindgen(getter)] pub fn stage(&self) -> Option<String> { match *self.inner { ck::HydrogenError::MissingSourceState { stage, .. } => Some(stage.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn state(&self) -> Option<String> { match *self.inner { ck::HydrogenError::MissingSourceState { state, .. } => Some(state.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::HydrogenError::ValenceAssignmentLength { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::ValenceAssignmentLength { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::ValenceAssignmentLength { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn current(&self) -> Option<u8> { match *self.inner { ck::HydrogenError::ExplicitHydrogenOverflow { current, .. } => Some(current), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::InvalidStereoTransition { bond, .. } => Some(bond.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn conformer(&self) -> Option<usize> { match *self.inner { ck::HydrogenError::CoordinatePlacement { conformer, .. } => Some(conformer), _ => None } }
    #[wasm_bindgen(getter)] pub fn dimension(&self) -> Option<String> { match *self.inner { ck::HydrogenError::CoordinatePlacement { dimension, .. } => Some(dimension.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<HydrogenTopology> { match self.inner.as_ref() { ck::HydrogenError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn coordinates(&self) -> Option<HydrogenCoordinates> { match self.inner.as_ref() { ck::HydrogenError::InvalidCoordinates(value) => Some(CoordinateValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = topologyEdit)] pub fn topology_edit(&self) -> Option<HydrogenEdit> { match self.inner.as_ref() { ck::HydrogenError::TopologyEdit(value) => Some(TopologyEditError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn mapping(&self) -> Option<HydrogenMapping> { match self.inner.as_ref() { ck::HydrogenError::InvalidMapping(value) => Some(MappingValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = queryState)] pub fn query_state(&self) -> Option<HydrogenQuery> { match self.inner.as_ref() { ck::HydrogenError::InvalidQueryState(value) => Some(QueryStateError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<HydrogenValence> { match self.inner.as_ref() { ck::HydrogenError::Valence(value) => Some(ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn sanitize(&self) -> Option<HydrogenSanitize> { match self.inner.as_ref() { ck::HydrogenError::Sanitize(value) => Some(SanitizeError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = HydrogenCause, cause = cause);
