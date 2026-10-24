//! Legacy stereo wrappers forward real sources, while retaining each wrapper value.
use crate::{
    atropisomer_errors, cip_rank_errors, double_stereo_errors,
    errors::{enum_error, exception},
    model_errors::{
        self, AtomPropertyError, BondValueError, StereoGroupError, TopologyValidationError,
    },
    potential_stereo_errors, ring_errors, stereo_order_errors, valence_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | ValenceError | ConjugationError | HybridizationError | RingFindingError | TopologyValidationError | CoordinateValidationError | BondValueError | StereoOrderError | AtomPropertyError | PropertyIntReadError | TemplateAttachmentOrderError | MappingValidationError | PropertyUIntReadError | CanonicalRankError | DoubleBondStereoError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type LegacyStereoCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type LegacyTopology;
    #[wasm_bindgen(typescript_type = "CipRankError")]
    pub type LegacyRank;
    #[wasm_bindgen(typescript_type = "DoubleBondStereoError")]
    pub type LegacyDouble;
    #[wasm_bindgen(typescript_type = "StereoOrderError")]
    pub type LegacyOrder;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type LegacyValence;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type LegacyAtomProperty;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type LegacyBondValue;
    #[wasm_bindgen(typescript_type = "PotentialStereoError")]
    pub type LegacyPotential;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type LegacyRings;
    #[wasm_bindgen(typescript_type = "AtropisomerError")]
    pub type LegacyAtropisomer;
}
pub(crate) fn cause(value: &ck::LegacyStereoError) -> Option<LegacyStereoCause> {
    match value {
        ck::LegacyStereoError::StereoGroup(value) => {
            Some(StereoGroupError::from_rust(*value).unchecked_into())
        }
        ck::LegacyStereoError::InvalidTopology(value) => {
            model_errors::topology_cause(value).map(JsCast::unchecked_into)
        }
        ck::LegacyStereoError::CipRank(value) => {
            cip_rank_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::LegacyStereoError::DoubleBond(value) => {
            double_stereo_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::LegacyStereoError::StereoOrder(value) => {
            stereo_order_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::LegacyStereoError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::LegacyStereoError::AtomProperty(_) | ck::LegacyStereoError::BondValue(_) => None,
        ck::LegacyStereoError::PotentialStereo(value) => {
            potential_stereo_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::LegacyStereoError::RingFinding(value) => {
            ring_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::LegacyStereoError::Atropisomer(value) => {
            atropisomer_errors::cause(value).map(JsCast::unchecked_into)
        }
    }
}
enum_error!(LegacyStereoError, LegacyStereoErrorKind, ck::LegacyStereoError, {
    StereoGroup => ck::LegacyStereoError::StereoGroup(_), InvalidTopology => ck::LegacyStereoError::InvalidTopology(_), CipRank => ck::LegacyStereoError::CipRank(_),
    DoubleBond => ck::LegacyStereoError::DoubleBond(_), StereoOrder => ck::LegacyStereoError::StereoOrder(_), Valence => ck::LegacyStereoError::Valence(_),
    AtomProperty => ck::LegacyStereoError::AtomProperty(_), BondValue => ck::LegacyStereoError::BondValue(_), PotentialStereo => ck::LegacyStereoError::PotentialStereo(_),
    RingFinding => ck::LegacyStereoError::RingFinding(_), Atropisomer => ck::LegacyStereoError::Atropisomer(_),
}, {
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<LegacyTopology> { match self.inner.as_ref() { ck::LegacyStereoError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = cipRank)] pub fn cip_rank(&self) -> Option<LegacyRank> { match self.inner.as_ref() { ck::LegacyStereoError::CipRank(value) => Some(cip_rank_errors::CipRankError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = doubleBond)] pub fn double_bond(&self) -> Option<LegacyDouble> { match self.inner.as_ref() { ck::LegacyStereoError::DoubleBond(value) => Some(double_stereo_errors::DoubleBondStereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = stereoOrder)] pub fn stereo_order(&self) -> Option<LegacyOrder> { match self.inner.as_ref() { ck::LegacyStereoError::StereoOrder(value) => Some(stereo_order_errors::StereoOrderError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<LegacyValence> { match self.inner.as_ref() { ck::LegacyStereoError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomProperty)] pub fn atom_property(&self) -> Option<LegacyAtomProperty> { match self.inner.as_ref() { ck::LegacyStereoError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = bondValue)] pub fn bond_value(&self) -> Option<LegacyBondValue> { match self.inner.as_ref() { ck::LegacyStereoError::BondValue(value) => Some(BondValueError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = potentialStereo)] pub fn potential_stereo(&self) -> Option<LegacyPotential> { match self.inner.as_ref() { ck::LegacyStereoError::PotentialStereo(value) => Some(potential_stereo_errors::PotentialStereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = ringFinding)] pub fn ring_finding(&self) -> Option<LegacyRings> { match self.inner.as_ref() { ck::LegacyStereoError::RingFinding(value) => Some(ring_errors::RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn atropisomer(&self) -> Option<LegacyAtropisomer> { match self.inner.as_ref() { ck::LegacyStereoError::Atropisomer(value) => Some(atropisomer_errors::AtropisomerError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = LegacyStereoCause, cause = cause);
