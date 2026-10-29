//! Scaffold diagnostics retain transparent native payloads and exact source edges.
use crate::{
    errors::{enum_error, exception},
    legacy_stereo_errors, matrix_errors,
    model_errors::{
        self, AtomPropertyError, BondValueError, MoleculePropertyError, TopologyValidationError,
    },
    property_errors, ring_errors, sanitize_leaf_errors, topology_errors, valence_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | ValenceError | ConjugationError | HybridizationError | RingFindingError | TopologyValidationError | CoordinateValidationError | BondValueError | StereoOrderError | AtomPropertyError | PropertyIntReadError | TemplateAttachmentOrderError | MappingValidationError | PropertyUIntReadError | CanonicalRankError | DoubleBondStereoError | MoleculePropertyError | UIntLexicalReadError"
    )]
    #[derive(Clone)]
    pub type ScaffoldCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type ScaffoldTopology;
    #[wasm_bindgen(typescript_type = "TopologyEditError")]
    pub type ScaffoldEdit;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type ScaffoldValence;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type ScaffoldRings;
    #[wasm_bindgen(typescript_type = "RadicalError")]
    pub type ScaffoldRadicals;
    #[wasm_bindgen(typescript_type = "LegacyStereoError")]
    pub type ScaffoldStereo;
    #[wasm_bindgen(typescript_type = "MoleculePropertyError")]
    pub type ScaffoldProperty;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type ScaffoldAtomProperty;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type ScaffoldBondProperty;
    #[wasm_bindgen(typescript_type = "PropertyUIntReadError")]
    pub type ScaffoldUInt;
    #[wasm_bindgen(typescript_type = "PropertyStringError")]
    pub type ScaffoldText;
    #[wasm_bindgen(typescript_type = "MatrixError")]
    pub type ScaffoldMatrix;
}
pub(crate) fn cause(value: &ck::ScaffoldError) -> Option<ScaffoldCause> {
    match value {
        ck::ScaffoldError::Topology(value) => {
            model_errors::topology_cause(value).map(JsCast::unchecked_into)
        }
        ck::ScaffoldError::Edit(value) => topology_errors::cause(value).map(JsCast::unchecked_into),
        ck::ScaffoldError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::ScaffoldError::Rings(value) => ring_errors::cause(value).map(JsCast::unchecked_into),
        ck::ScaffoldError::Radicals(value) => {
            sanitize_leaf_errors::radical_cause(value).map(JsCast::unchecked_into)
        }
        ck::ScaffoldError::Stereo(value) => {
            legacy_stereo_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::ScaffoldError::UIntProperty(ck::PropertyUIntReadError::Lexical { source, .. }) => {
            Some(property_errors::UIntLexicalReadError::from_rust(source.clone()).unchecked_into())
        }
        ck::ScaffoldError::Matrix(value) => matrix_errors::cause(value).map(JsCast::unchecked_into),
        ck::ScaffoldError::Property(_)
        | ck::ScaffoldError::AtomProperty(_)
        | ck::ScaffoldError::BondProperty(_)
        | ck::ScaffoldError::TextProperty(_)
        | ck::ScaffoldError::UIntProperty(_)
        | ck::ScaffoldError::MissingImplicitHydrogens { .. }
        | ck::ScaffoldError::ValenceRowCount { .. }
        | ck::ScaffoldError::MissingRingInfo
        | ck::ScaffoldError::MatrixSize { .. } => None,
    }
}
enum_error!(ScaffoldError, ScaffoldErrorKind, ck::ScaffoldError, {
    Topology => ck::ScaffoldError::Topology(_), Edit => ck::ScaffoldError::Edit(_), Valence => ck::ScaffoldError::Valence(_), Rings => ck::ScaffoldError::Rings(_),
    Radicals => ck::ScaffoldError::Radicals(_), Stereo => ck::ScaffoldError::Stereo(_), Property => ck::ScaffoldError::Property(_), AtomProperty => ck::ScaffoldError::AtomProperty(_),
    BondProperty => ck::ScaffoldError::BondProperty(_), UIntProperty => ck::ScaffoldError::UIntProperty(_), TextProperty => ck::ScaffoldError::TextProperty(_), Matrix => ck::ScaffoldError::Matrix(_),
    MissingImplicitHydrogens => ck::ScaffoldError::MissingImplicitHydrogens { .. }, ValenceRowCount => ck::ScaffoldError::ValenceRowCount { .. }, MissingRingInfo => ck::ScaffoldError::MissingRingInfo, MatrixSize => ck::ScaffoldError::MatrixSize { .. },
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::ScaffoldError::MissingImplicitHydrogens { atom } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::ScaffoldError::ValenceRowCount { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::ScaffoldError::ValenceRowCount { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::ScaffoldError::ValenceRowCount { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn atoms(&self) -> Option<usize> { match *self.inner { ck::ScaffoldError::MatrixSize { atoms } => Some(atoms), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<ScaffoldTopology> { match self.inner.as_ref() { ck::ScaffoldError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn edit(&self) -> Option<ScaffoldEdit> { match self.inner.as_ref() { ck::ScaffoldError::Edit(value) => Some(topology_errors::TopologyEditError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<ScaffoldValence> { match self.inner.as_ref() { ck::ScaffoldError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn rings(&self) -> Option<ScaffoldRings> { match self.inner.as_ref() { ck::ScaffoldError::Rings(value) => Some(ring_errors::RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn radicals(&self) -> Option<ScaffoldRadicals> { match self.inner.as_ref() { ck::ScaffoldError::Radicals(value) => Some(sanitize_leaf_errors::RadicalError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn stereo(&self) -> Option<ScaffoldStereo> { match self.inner.as_ref() { ck::ScaffoldError::Stereo(value) => Some(legacy_stereo_errors::LegacyStereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = propertyError)] pub fn property_error(&self) -> Option<ScaffoldProperty> { match self.inner.as_ref() { ck::ScaffoldError::Property(value) => Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomProperty)] pub fn atom_property(&self) -> Option<ScaffoldAtomProperty> { match self.inner.as_ref() { ck::ScaffoldError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = bondProperty)] pub fn bond_property(&self) -> Option<ScaffoldBondProperty> { match self.inner.as_ref() { ck::ScaffoldError::BondProperty(value) => Some(BondValueError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = uintProperty)] pub fn uint_property(&self) -> Option<ScaffoldUInt> { match self.inner.as_ref() { ck::ScaffoldError::UIntProperty(value) => Some(property_errors::PropertyUIntReadError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = textProperty)] pub fn text_property(&self) -> Option<ScaffoldText> { match self.inner.as_ref() { ck::ScaffoldError::TextProperty(value) => Some(property_errors::PropertyStringError::from_rust(*value).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn matrix(&self) -> Option<ScaffoldMatrix> { match self.inner.as_ref() { ck::ScaffoldError::Matrix(value) => Some(matrix_errors::MatrixError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = ScaffoldCause, cause = cause);
