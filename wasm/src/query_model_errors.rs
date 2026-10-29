//! Query-model validation payloads retain the actual source annotations.
use crate::{
    errors::{enum_error, exception},
    model_errors::*,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "StereoGroupError | TemplateAttachmentOrderError")]
    #[derive(Clone)]
    pub type QueryGraphCause;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type QueryCoordinatePayload;
}
pub(crate) fn cause(e: &ck::QueryGraphError) -> Option<QueryGraphCause> {
    match e {
        ck::QueryGraphError::StereoGroup(v) => {
            Some(StereoGroupError::from_rust(*v).unchecked_into())
        }
        ck::QueryGraphError::TemplateAttachmentOrder { source, .. } => {
            Some(TemplateAttachmentOrderError::from_rust(source.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(QueryGraphError,QueryGraphErrorKind,ck::QueryGraphError,{
AtomIdMismatch=>ck::QueryGraphError::AtomIdMismatch{..},
TemplateAttachmentOrder=>ck::QueryGraphError::TemplateAttachmentOrder{..},
BondIdMismatch=>ck::QueryGraphError::BondIdMismatch{..},
StereoAtomOutOfRange=>ck::QueryGraphError::StereoAtomOutOfRange{..},
StereoGroupAtomOutOfRange=>ck::QueryGraphError::StereoGroupAtomOutOfRange{..},
StereoGroupBondOutOfRange=>ck::QueryGraphError::StereoGroupBondOutOfRange{..},
SubstanceGroupIdMismatch=>ck::QueryGraphError::SubstanceGroupIdMismatch{..},
SubstanceGroupAtomOutOfRange=>ck::QueryGraphError::SubstanceGroupAtomOutOfRange{..},
SubstanceGroupBondOutOfRange=>ck::QueryGraphError::SubstanceGroupBondOutOfRange{..},
SubstanceGroupParentOutOfRange=>ck::QueryGraphError::SubstanceGroupParentOutOfRange{..},
StereoGroup=>ck::QueryGraphError::StereoGroup(_),InvalidBondEndpoint=>ck::QueryGraphError::InvalidBondEndpoint(_),CoordinateValidation=>ck::QueryGraphError::CoordinateValidation(_),AdjacencyMismatch=>ck::QueryGraphError::AdjacencyMismatch,
},{
#[wasm_bindgen(getter,js_name=position)] pub fn position(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::AtomIdMismatch{position,..}=>Some(*position),
ck::QueryGraphError::BondIdMismatch{position,..}=>Some(*position),
ck::QueryGraphError::SubstanceGroupIdMismatch{position,..}=>Some(*position),
_=>None}}
#[wasm_bindgen(getter,js_name=id)] pub fn id(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::AtomIdMismatch{id,..}=>Some(id.index()),
ck::QueryGraphError::BondIdMismatch{id,..}=>Some(id.index()),
ck::QueryGraphError::SubstanceGroupIdMismatch{id,..}=>Some(id.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=atom)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::TemplateAttachmentOrder{atom,..}=>Some(atom.index()),
ck::QueryGraphError::StereoGroupAtomOutOfRange{atom,..}=>Some(atom.index()),
ck::QueryGraphError::SubstanceGroupAtomOutOfRange{atom,..}=>Some(atom.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=bond)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::StereoAtomOutOfRange{bond,..}=>Some(bond.index()),
ck::QueryGraphError::StereoGroupBondOutOfRange{bond,..}=>Some(bond.index()),
ck::QueryGraphError::SubstanceGroupBondOutOfRange{bond,..}=>Some(bond.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=begin)] pub fn begin(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::StereoAtomOutOfRange{begin,..}=>Some(begin.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=end)] pub fn end(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::StereoAtomOutOfRange{end,..}=>Some(end.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::StereoAtomOutOfRange{atom_count,..}=>Some(*atom_count),
ck::QueryGraphError::StereoGroupAtomOutOfRange{atom_count,..}=>Some(*atom_count),
ck::QueryGraphError::SubstanceGroupAtomOutOfRange{atom_count,..}=>Some(*atom_count),
_=>None}}
#[wasm_bindgen(getter,js_name=bondCount)] pub fn bond_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::StereoGroupBondOutOfRange{bond_count,..}=>Some(*bond_count),
ck::QueryGraphError::SubstanceGroupBondOutOfRange{bond_count,..}=>Some(*bond_count),
_=>None}}
#[wasm_bindgen(getter,js_name=sgroup)] pub fn sgroup(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::SubstanceGroupAtomOutOfRange{sgroup,..}=>Some(sgroup.index()),
ck::QueryGraphError::SubstanceGroupBondOutOfRange{sgroup,..}=>Some(sgroup.index()),
ck::QueryGraphError::SubstanceGroupParentOutOfRange{sgroup,..}=>Some(sgroup.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=parent)] pub fn parent(&self)->Option<usize>{match self.inner.as_ref(){
ck::QueryGraphError::SubstanceGroupParentOutOfRange{parent,..}=>Some(parent.index()),
_=>None}}
#[wasm_bindgen(getter,js_name=invalidBondEndpoint)] pub fn invalid_bond_endpoint(&self)->Option<usize>{match self.inner.as_ref(){ck::QueryGraphError::InvalidBondEndpoint(v)=>Some(*v),_=>None}}
#[wasm_bindgen(getter)] pub fn source(&self)->Option<QueryGraphCause>{cause(self.inner.as_ref())}
#[wasm_bindgen(getter,js_name=coordinateValidation)] pub fn coordinate_validation(&self)->Option<QueryCoordinatePayload>{match self.inner.as_ref(){ck::QueryGraphError::CoordinateValidation(v)=>Some(CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
},source=QueryGraphCause,cause=cause);
enum_error!(QueryAtomConversionError,QueryAtomConversionErrorKind,ck::QueryAtomConversionError,{NonElementAtomicNumber=>ck::QueryAtomConversionError::NonElementAtomicNumber{..}},{
#[wasm_bindgen(getter)] pub fn atom(&self)->usize {let ck::QueryAtomConversionError::NonElementAtomicNumber{atom,..}=self.inner.as_ref();atom.index()}
#[wasm_bindgen(getter,js_name=atomicNumber)] pub fn atomic_number(&self)->u8 {let ck::QueryAtomConversionError::NonElementAtomicNumber{atomic_number,..}=self.inner.as_ref();*atomic_number}
});
