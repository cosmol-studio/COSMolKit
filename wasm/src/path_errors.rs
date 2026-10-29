//! Native path diagnostics shared by graph queries and depiction.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | MoleculePropertyError | BondValueError | AtomPropertyError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type PathCause;
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | MoleculePropertyError | BondValueError | AtomPropertyError | TopologyValidationError | PeriodicTableError | QueryGraphError | MappingValidationError | MatrixError | TemplateAttachmentOrderError"
    )]
    pub type PathPayload;
}
pub(crate) fn source(e: &ck::PathError) -> Option<PathCause> {
    use ck::PathError as E;
    match e {
        E::StereoGroup(v) => {
            Some(crate::model_errors::StereoGroupError::from_rust(*v).unchecked_into())
        }
        E::MoleculeProperty(v) => {
            Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into())
        }
        E::BondProperty(v) => {
            Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into())
        }
        E::AtomProperty(v) => {
            Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into())
        }
        E::InvalidTopology(v) => crate::model_errors::topology_cause(v).map(JsCast::unchecked_into),
        E::QueryGraph(v) => crate::query_model_errors::cause(v).map(JsCast::unchecked_into),
        E::Matrix(v) => crate::matrix_errors::cause(v).map(JsCast::unchecked_into),
        E::TemplateAttachmentRemap { source, .. } => Some(
            crate::model_errors::TemplateAttachmentOrderError::from_rust(source.clone())
                .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(PathError,PathErrorKind,ck::PathError,{
    StereoGroup=>ck::PathError::StereoGroup(_),MoleculeProperty=>ck::PathError::MoleculeProperty(_),BondProperty=>ck::PathError::BondProperty(_),AtomProperty=>ck::PathError::AtomProperty(_),InvalidTopology=>ck::PathError::InvalidTopology(_),
    AtomOutOfRange=>ck::PathError::AtomOutOfRange{..},EqualShortestPathEndpoints=>ck::PathError::EqualShortestPathEndpoints{..},InvalidLengthRange=>ck::PathError::InvalidLengthRange{..},LengthOverflow=>ck::PathError::LengthOverflow{..},AtomPathOutOfRange=>ck::PathError::AtomPathOutOfRange{..},IgnoredAtomMaskLength=>ck::PathError::IgnoredAtomMaskLength{..},ExtraInvariantLength=>ck::PathError::ExtraInvariantLength{..},PeriodicTable=>ck::PathError::PeriodicTable(_),InvalidSubsetTopology=>ck::PathError::InvalidSubsetTopology(_),QueryGraph=>ck::PathError::QueryGraph(_),Mapping=>ck::PathError::Mapping(_),TemplateAttachmentRemap=>ck::PathError::TemplateAttachmentRemap{..},Matrix=>ck::PathError::Matrix(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<PathPayload>{use ck::PathError as E;Some(match self.inner.as_ref(){
        E::StereoGroup(v)=>crate::model_errors::StereoGroupError::from_rust(*v),E::MoleculeProperty(v)=>crate::model_errors::MoleculePropertyError::from_rust(v.clone()),E::BondProperty(v)=>crate::model_errors::BondValueError::from_rust(v.clone()),E::AtomProperty(v)=>crate::model_errors::AtomPropertyError::from_rust(v.clone()),
        E::InvalidTopology(v)|E::InvalidSubsetTopology(v)=>crate::model_errors::TopologyValidationError::from_rust(v.clone()),E::PeriodicTable(v)=>crate::periodic_errors::PeriodicTableError::from_rust(v.clone()),E::QueryGraph(v)=>crate::query_model_errors::QueryGraphError::from_rust(v.clone()),E::Mapping(v)=>crate::model_errors::MappingValidationError::from_rust(v.clone()),E::Matrix(v)=>crate::matrix_errors::MatrixError::from_rust(v.clone()),E::TemplateAttachmentRemap{source,..}=>crate::model_errors::TemplateAttachmentOrderError::from_rust(source.clone()),_=>return None,
    }.unchecked_into())}
    #[wasm_bindgen(getter)] pub fn role(&self)->Option<String>{match self.inner.as_ref(){ck::PathError::AtomOutOfRange{role,..}=>Some((*role).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{use ck::PathError as E;match self.inner.as_ref(){E::AtomOutOfRange{atom,..}|E::EqualShortestPathEndpoints{atom}|E::AtomPathOutOfRange{atom,..}=>Some(atom.index()),_=>None}}
    #[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{use ck::PathError as E;match self.inner.as_ref(){E::AtomOutOfRange{atom_count,..}|E::AtomPathOutOfRange{atom_count,..}=>Some(*atom_count),_=>None}}
    #[wasm_bindgen(getter)] pub fn lower(&self)->Option<usize>{match self.inner.as_ref(){ck::PathError::InvalidLengthRange{lower,..}=>Some(*lower),_=>None}}
    #[wasm_bindgen(getter)] pub fn upper(&self)->Option<usize>{match self.inner.as_ref(){ck::PathError::InvalidLengthRange{upper,..}=>Some(*upper),_=>None}}
    #[wasm_bindgen(getter)] pub fn length(&self)->Option<usize>{match self.inner.as_ref(){ck::PathError::LengthOverflow{length}=>Some(*length),_=>None}}
    #[wasm_bindgen(getter)] pub fn position(&self)->Option<usize>{match self.inner.as_ref(){ck::PathError::AtomPathOutOfRange{position,..}=>Some(*position),_=>None}}
    #[wasm_bindgen(getter)] pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::PathError::IgnoredAtomMaskLength{expected,..}|ck::PathError::ExtraInvariantLength{expected,..}=>Some(*expected),_=>None}}
    #[wasm_bindgen(getter)] pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::PathError::IgnoredAtomMaskLength{actual,..}|ck::PathError::ExtraInvariantLength{actual,..}=>Some(*actual),_=>None}}
    #[wasm_bindgen(getter)] pub fn carrier(&self)->Option<usize>{match self.inner.as_ref(){ck::PathError::TemplateAttachmentRemap{carrier,..}=>Some(carrier.index()),_=>None}}
},source=PathCause,cause=source);
