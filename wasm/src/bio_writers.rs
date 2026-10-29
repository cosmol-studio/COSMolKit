//! BIO text writers preserve the canonical error graph and defaults.
use crate::{
    bio_hierarchy::BioStructure,
    bio_read_errors::{BioStructureError, CifReadError},
    errors::{enum_error, exception},
    external_errors::IoError,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
ck::__ck_configuration_BioPdbWriteParams!(crate::project_parameters, ck::BioPdbWriteParams, {}, {
    BioPdbWriteParams
});
ck::__ck_configuration_BioMmcifWriteParams!(
    crate::project_parameters,
    ck::BioMmcifWriteParams,
    {},
    { BioMmcifWriteParams }
);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "BioStructureError | IoError")]
    #[derive(Clone)]
    pub type BioPdbWriteCause;
    #[wasm_bindgen(typescript_type = "BioStructureError | CifReadError | IoError")]
    #[derive(Clone)]
    pub type BioMmcifWriteCause;
}
enum_error!(BioPdbWriteError,BioPdbWriteErrorKind,ck::BioPdbWriteError,{InvalidStructure=>ck::BioPdbWriteError::InvalidStructure(_),ChainNameTooLong=>ck::BioPdbWriteError::ChainNameTooLong{..},NegativeSerial=>ck::BioPdbWriteError::NegativeSerial{..},SerialIncrementOverflow=>ck::BioPdbWriteError::SerialIncrementOverflow{..},SerialOffsetOverflow=>ck::BioPdbWriteError::SerialOffsetOverflow{..},NegativeBase36Value=>ck::BioPdbWriteError::NegativeBase36Value{..},UnrepresentableField=>ck::BioPdbWriteError::UnrepresentableField{..},Io=>ck::BioPdbWriteError::Io{..}},{
#[wasm_bindgen(getter)]pub fn chain(&self)->Option<String>{match self.inner.as_ref(){ck::BioPdbWriteError::ChainNameTooLong{chain,..}=>Some(chain.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn length(&self)->Option<usize>{match self.inner.as_ref(){ck::BioPdbWriteError::ChainNameTooLong{length,..}=>Some(*length),_=>None}}
#[wasm_bindgen(getter)]pub fn serial(&self)->Option<i32>{match self.inner.as_ref(){ck::BioPdbWriteError::NegativeSerial{serial,..}|ck::BioPdbWriteError::SerialOffsetOverflow{serial,..}=>Some(*serial),_=>None}}
#[wasm_bindgen(getter)]pub fn from(&self)->Option<i32>{match self.inner.as_ref(){ck::BioPdbWriteError::SerialIncrementOverflow{from,..}=>Some(*from),_=>None}}
#[wasm_bindgen(getter)]pub fn value(&self)->Option<i32>{match self.inner.as_ref(){ck::BioPdbWriteError::NegativeBase36Value{value,..}=>Some(*value),_=>None}}
#[wasm_bindgen(getter)]pub fn field(&self)->Option<String>{match self.inner.as_ref(){ck::BioPdbWriteError::UnrepresentableField{field,..}=>Some(field.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn detail(&self)->Option<String>{match self.inner.as_ref(){ck::BioPdbWriteError::UnrepresentableField{detail,..}=>Some(detail.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn path(&self)->Option<String>{match self.inner.as_ref(){ck::BioPdbWriteError::Io{path,..}=>Some(path.display().to_string()),_=>None}}
},source=BioPdbWriteCause,cause=|owner:&ReadonlyValue<ck::BioPdbWriteError>|Some(match owner.as_ref(){ck::BioPdbWriteError::InvalidStructure(_)=>BioStructureError::from_owner(owner.child(|v|match v{ck::BioPdbWriteError::InvalidStructure(v)=>v,_=>unreachable!("immutable structure source")})),ck::BioPdbWriteError::Io{..}=>IoError::from_owner(owner.child(|v|match v{ck::BioPdbWriteError::Io{source,..}=>source,_=>unreachable!("immutable IO source")})),_=>return None}.unchecked_into()));
enum_error!(BioMmcifWriteError,BioMmcifWriteErrorKind,ck::BioMmcifWriteError,{Structure=>ck::BioMmcifWriteError::Structure(_),Cif=>ck::BioMmcifWriteError::Cif(_),InvalidText=>ck::BioMmcifWriteError::InvalidText{..},Io=>ck::BioMmcifWriteError::Io(_),FileWrite=>ck::BioMmcifWriteError::FileWrite{..}},{
#[wasm_bindgen(getter)]pub fn field(&self)->Option<String>{match self.inner.as_ref(){ck::BioMmcifWriteError::InvalidText{field}=>Some(field.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn path(&self)->Option<String>{match self.inner.as_ref(){ck::BioMmcifWriteError::FileWrite{path,..}=>Some(path.display().to_string()),_=>None}}
},source=BioMmcifWriteCause,cause=|owner:&ReadonlyValue<ck::BioMmcifWriteError>|Some(match owner.as_ref(){
ck::BioMmcifWriteError::Structure(_)=>BioStructureError::from_owner(owner.child(|v|match v{ck::BioMmcifWriteError::Structure(v)=>v,_=>unreachable!("immutable structure source")})),
ck::BioMmcifWriteError::Cif(_)=>CifReadError::from_owner(owner.child(|v|match v{ck::BioMmcifWriteError::Cif(v)=>v,_=>unreachable!("immutable CIF source")})),
ck::BioMmcifWriteError::Io(_)|ck::BioMmcifWriteError::FileWrite{..}=>IoError::from_owner(owner.child(|v|match v{ck::BioMmcifWriteError::Io(v)|ck::BioMmcifWriteError::FileWrite{source:v,..}=>v,_=>unreachable!("immutable IO source")})),
ck::BioMmcifWriteError::InvalidText{..}=>return None}.unchecked_into()));
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=toPdbWithParams)]
    pub fn to_pdb_with_params(&self, params: BioPdbWriteParamsInput) -> Result<String, JsValue> {
        self.to_pdb(Some(params))
    }
    #[wasm_bindgen(js_name=toMmcifWithParams)]
    pub fn to_mmcif_with_params(
        &self,
        params: BioMmcifWriteParamsInput,
    ) -> Result<String, JsValue> {
        self.to_mmcif(Some(params))
    }
    #[wasm_bindgen(js_name=toPdb)]
    pub fn to_pdb(&self, params: Option<BioPdbWriteParamsInput>) -> Result<String, JsValue> {
        BioPdbWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_pdb_with_params(p)
                .map_err(BioPdbWriteError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toMmcif)]
    pub fn to_mmcif(&self, params: Option<BioMmcifWriteParamsInput>) -> Result<String, JsValue> {
        BioMmcifWriteParams::with_input(params.as_ref(), |p| {
            self.0
                .to_mmcif_with_params(p)
                .map_err(BioMmcifWriteError::from_rust)
        })
    }
}
