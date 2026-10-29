//! IO finalization dependencies are exact native error projections.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | CoordinateValidationError | QueryStateError | AtomPropertyError | TopologyEditError | ValenceError | HydrogenError"
    )]
    #[derive(Clone)]
    pub type AttachmentCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type AttachmentTopology;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type AttachmentCoordinates;
    #[wasm_bindgen(typescript_type = "QueryStateError")]
    pub type AttachmentQuery;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type AttachmentProperty;
    #[wasm_bindgen(typescript_type = "TopologyEditError")]
    pub type AttachmentEdit;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type AttachmentValence;
    #[wasm_bindgen(typescript_type = "HydrogenError")]
    pub type AttachmentPlacement;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type DirectionValence;
    #[wasm_bindgen(typescript_type = "AtropisomerError")]
    pub type DirectionAtrop;
}
fn attachment_cause(e: &ck::AttachmentExpansionError) -> Option<AttachmentCause> {
    match e {
        ck::AttachmentExpansionError::Topology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::AttachmentExpansionError::Coordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::AttachmentExpansionError::Query(v) => {
            Some(crate::query_errors::QueryStateError::from_rust(v.clone()).unchecked_into())
        }
        ck::AttachmentExpansionError::Property(v) => {
            Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::AttachmentExpansionError::Edit(v) => {
            Some(crate::topology_errors::TopologyEditError::from_rust(v.clone()).unchecked_into())
        }
        ck::AttachmentExpansionError::Valence(v) => {
            Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into())
        }
        ck::AttachmentExpansionError::Placement(v) => {
            Some(crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into())
        }
        ck::AttachmentExpansionError::ValueCount { .. } => None,
    }
}
enum_error!(AttachmentExpansionError,AttachmentExpansionErrorKind,ck::AttachmentExpansionError,{
Topology=>ck::AttachmentExpansionError::Topology(_),
Coordinates=>ck::AttachmentExpansionError::Coordinates(_),
Query=>ck::AttachmentExpansionError::Query(_),
Property=>ck::AttachmentExpansionError::Property(_),
Edit=>ck::AttachmentExpansionError::Edit(_),
Valence=>ck::AttachmentExpansionError::Valence(_),
Placement=>ck::AttachmentExpansionError::Placement(_),
ValueCount=>ck::AttachmentExpansionError::ValueCount{..},},{
#[wasm_bindgen(getter)] pub fn topology(&self)->Option<AttachmentTopology>{match self.inner.as_ref(){ck::AttachmentExpansionError::Topology(v)=>Some(crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn coordinates(&self)->Option<AttachmentCoordinates>{match self.inner.as_ref(){ck::AttachmentExpansionError::Coordinates(v)=>Some(crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn query(&self)->Option<AttachmentQuery>{match self.inner.as_ref(){ck::AttachmentExpansionError::Query(v)=>Some(crate::query_errors::QueryStateError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn property(&self)->Option<AttachmentProperty>{match self.inner.as_ref(){ck::AttachmentExpansionError::Property(v)=>Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn edit(&self)->Option<AttachmentEdit>{match self.inner.as_ref(){ck::AttachmentExpansionError::Edit(v)=>Some(crate::topology_errors::TopologyEditError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn valence(&self)->Option<AttachmentValence>{match self.inner.as_ref(){ck::AttachmentExpansionError::Valence(v)=>Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn placement(&self)->Option<AttachmentPlacement>{match self.inner.as_ref(){ck::AttachmentExpansionError::Placement(v)=>Some(crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::AttachmentExpansionError::ValueCount{actual,..}=>Some(*actual),_=>None}}
#[wasm_bindgen(getter)] pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::AttachmentExpansionError::ValueCount{expected,..}=>Some(*expected),_=>None}}
},source=AttachmentCause,cause=attachment_cause);
pub(crate) fn direction_cause(
    e: &ck::BondDirectionStereoError,
) -> Option<crate::atropisomer_errors::AtropisomerCause> {
    match e {
        ck::BondDirectionStereoError::Valence(v) => {
            crate::valence_errors::cause(v).map(JsCast::unchecked_into)
        }
        ck::BondDirectionStereoError::Atropisomer(v) => crate::atropisomer_errors::cause(v),
        ck::BondDirectionStereoError::InvalidState(_) => None,
    }
}
enum_error!(BondDirectionStereoError,BondDirectionStereoErrorKind,ck::BondDirectionStereoError,{InvalidState=>ck::BondDirectionStereoError::InvalidState(_),Valence=>ck::BondDirectionStereoError::Valence(_),Atropisomer=>ck::BondDirectionStereoError::Atropisomer(_)},{
#[wasm_bindgen(getter,js_name=invalidState)] pub fn invalid_state(&self)->Option<String>{match self.inner.as_ref(){ck::BondDirectionStereoError::InvalidState(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter)] pub fn valence(&self)->Option<DirectionValence>{match self.inner.as_ref(){ck::BondDirectionStereoError::Valence(v)=>Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)] pub fn atropisomer(&self)->Option<DirectionAtrop>{match self.inner.as_ref(){ck::BondDirectionStereoError::Atropisomer(v)=>Some(crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()).unchecked_into()),_=>None}}
},source=crate::atropisomer_errors::AtropisomerCause,cause=direction_cause);
