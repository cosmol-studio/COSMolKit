//! Wedge errors retain source-defined transparent cause edges.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | CoordinateValidationError | StereoGroupError | ValenceError | ConjugationError | HybridizationError | RingFindingError | BondValueError | StereoOrderError | AtomPropertyError | PropertyIntReadError | PropertyUIntReadError | CanonicalRankError | DoubleBondStereoError | TemplateAttachmentOrderError | RingDecomposerError"
    )]
    #[derive(Clone)]
    pub type WedgeCause;
    #[wasm_bindgen(typescript_type = "number | string")]
    pub type WedgeValue;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type WedgeTopology;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type WedgeCoordinates;
    #[wasm_bindgen(typescript_type = "StereoError")]
    pub type WedgeGeometry;
    #[wasm_bindgen(typescript_type = "StereoOrderError")]
    pub type WedgeStereoOrder;
    #[wasm_bindgen(typescript_type = "PotentialStereoError")]
    pub type WedgePotentialStereo;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type WedgeRingFinding;
    #[wasm_bindgen(typescript_type = "AtropisomerError")]
    pub type WedgeAtropisomer;
}
pub(crate) fn cause(e: &ck::WedgeError) -> Option<WedgeCause> {
    match e {
        ck::WedgeError::Topology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::WedgeError::Coordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::WedgeError::Geometry(v) => crate::stereo_errors::cause(v).map(JsCast::unchecked_into),
        ck::WedgeError::StereoOrder(v) => {
            crate::stereo_order_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::WedgeError::PotentialStereo(v) => {
            crate::potential_stereo_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::WedgeError::RingFinding(v) => crate::ring_errors::cause(v).map(JsCast::unchecked_into),
        ck::WedgeError::Atropisomer(v) => {
            crate::atropisomer_errors::cause(v).map(JsCast::unchecked_into)
        }
        _ => None,
    }
}
enum_error!(WedgeError,WedgeErrorKind,ck::WedgeError,{
Topology=>ck::WedgeError::Topology(_),
Coordinates=>ck::WedgeError::Coordinates(_),
Geometry=>ck::WedgeError::Geometry(_),
StereoOrder=>ck::WedgeError::StereoOrder(_),
PotentialStereo=>ck::WedgeError::PotentialStereo(_),
RingFinding=>ck::WedgeError::RingFinding(_),
Atropisomer=>ck::WedgeError::Atropisomer(_),
SourcePrecondition=>ck::WedgeError::SourcePrecondition{..},
SourceStateIndex=>ck::WedgeError::SourceStateIndex{..},
SourceIndexWidth=>ck::WedgeError::SourceIndexWidth{..},
SourceSignedScoreOverflow=>ck::WedgeError::SourceSignedScoreOverflow{..},
UnsignedRankOverflow=>ck::WedgeError::UnsignedRankOverflow{..},
BondOutOfRange=>ck::WedgeError::BondOutOfRange{..},
NonSingleBond=>ck::WedgeError::NonSingleBond{..},
CenterNotIncident=>ck::WedgeError::CenterNotIncident{..},
UnsupportedCenterTag=>ck::WedgeError::UnsupportedCenterTag{..},
InvalidStereoRank=>ck::WedgeError::InvalidStereoRank{..},
InvalidStereoRankType=>ck::WedgeError::InvalidStereoRankType{..},
},{
#[wasm_bindgen(getter,js_name=topology)] pub fn topology(&self)->Option<WedgeTopology>{match self.inner.as_ref(){ck::WedgeError::Topology(v)=>Some(crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=coordinates)] pub fn coordinates(&self)->Option<WedgeCoordinates>{match self.inner.as_ref(){ck::WedgeError::Coordinates(v)=>Some(crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=geometry)] pub fn geometry(&self)->Option<WedgeGeometry>{match self.inner.as_ref(){ck::WedgeError::Geometry(v)=>Some(crate::stereo_errors::StereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stereoOrder)] pub fn stereo_order(&self)->Option<WedgeStereoOrder>{match self.inner.as_ref(){ck::WedgeError::StereoOrder(v)=>Some(crate::stereo_order_errors::StereoOrderError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=potentialStereo)] pub fn potential_stereo(&self)->Option<WedgePotentialStereo>{match self.inner.as_ref(){ck::WedgeError::PotentialStereo(v)=>Some(crate::potential_stereo_errors::PotentialStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=ringFinding)] pub fn ring_finding(&self)->Option<WedgeRingFinding>{match self.inner.as_ref(){ck::WedgeError::RingFinding(v)=>Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atropisomer)] pub fn atropisomer(&self)->Option<WedgeAtropisomer>{match self.inner.as_ref(){ck::WedgeError::Atropisomer(v)=>Some(crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=message)] pub fn message(&self)->Option<String>{match self.inner.as_ref(){
ck::WedgeError::SourcePrecondition{message,..}=>Some(message.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=state)] pub fn state(&self)->Option<String>{match self.inner.as_ref(){
ck::WedgeError::SourceStateIndex{state,..}=>Some(state.to_string()),
ck::WedgeError::SourceIndexWidth{state,..}=>Some(state.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=index)] pub fn index(&self)->Option<usize>{match self.inner.as_ref(){
ck::WedgeError::SourceStateIndex{index,..}=>Some(*index),
ck::WedgeError::SourceIndexWidth{index,..}=>Some(*index),
_=>None}}
#[wasm_bindgen(getter,js_name=count)] pub fn count(&self)->Option<usize>{match self.inner.as_ref(){
ck::WedgeError::SourceStateIndex{count,..}=>Some(*count),
_=>None}}
#[wasm_bindgen(getter,js_name=atom)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){
ck::WedgeError::SourceSignedScoreOverflow{atom,..}=>Some(atom.index()),
ck::WedgeError::UnsignedRankOverflow{atom,..}=>Some(atom.index()),
ck::WedgeError::InvalidStereoRank{atom,..}=>Some(atom.index()),
ck::WedgeError::InvalidStereoRankType{atom,..}=>Some(atom.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=operation)] pub fn operation(&self)->Option<String>{match self.inner.as_ref(){
ck::WedgeError::SourceSignedScoreOverflow{operation,..}=>Some(operation.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=property)] pub fn property(&self)->Option<String>{match self.inner.as_ref(){
ck::WedgeError::UnsignedRankOverflow{property,..}=>Some(property.to_string()),
ck::WedgeError::InvalidStereoRank{property,..}=>Some(property.to_string()),
ck::WedgeError::InvalidStereoRankType{property,..}=>Some(property.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=bond)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){
ck::WedgeError::BondOutOfRange{bond,..}=>Some(bond.index()),
ck::WedgeError::NonSingleBond{bond,..}=>Some(bond.index()),
ck::WedgeError::CenterNotIncident{bond,..}=>Some(bond.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=bondCount)] pub fn bond_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::WedgeError::BondOutOfRange{bond_count,..}=>Some(*bond_count),
_=>None}}
#[wasm_bindgen(getter,js_name=center)] pub fn center(&self)->Option<usize>{match self.inner.as_ref(){
ck::WedgeError::CenterNotIncident{center,..}=>Some(center.index()),
ck::WedgeError::UnsupportedCenterTag{center,..}=>Some(center.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=tag)] pub fn tag(&self)->Option<crate::enums::ChiralTag>{match self.inner.as_ref(){
ck::WedgeError::UnsupportedCenterTag{tag,..}=>Some((*tag).into()),
_=>None}}
#[wasm_bindgen(getter,js_name=kind)] pub fn kind(&self)->Option<crate::model_errors::PropertyValueKind>{match self.inner.as_ref(){
ck::WedgeError::InvalidStereoRankType{kind,..}=>Some((*kind).into()),
_=>None}}
#[wasm_bindgen(getter)] pub fn value(&self)->Result<Option<WedgeValue>,JsValue>{match self.inner.as_ref(){ck::WedgeError::UnsignedRankOverflow{value,..}=>Ok(Some(JsValue::from(*value).unchecked_into())),ck::WedgeError::InvalidStereoRank{value,..}=>crate::property_errors::text(value).map(|v|Some(JsValue::from(v).unchecked_into())),_=>Ok(None)}}
},source=WedgeCause,cause=cause);
