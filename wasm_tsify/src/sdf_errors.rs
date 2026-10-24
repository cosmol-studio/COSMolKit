//! Structured SDF read failures, including record framing offsets.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | MolPostError | TopologyValidationError | CoordinateValidationError | QueryGraphError | AtomPropertyError | BondValueError | MoleculePropertyError | SdfReadError"
    )]
    #[derive(Clone)]
    pub type SdfReadErrorCause;
    #[wasm_bindgen(typescript_type = "StereoGroupError")]
    pub type SdfReadErrorStereoGroup;
    #[wasm_bindgen(typescript_type = "MolPostError")]
    pub type SdfReadErrorMolPost;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type SdfReadErrorTopology;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type SdfReadErrorCoordinates;
    #[wasm_bindgen(typescript_type = "QueryGraphError")]
    pub type SdfReadErrorQueryGraph;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type SdfReadErrorAtomProperty;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type SdfReadErrorBondValue;
    #[wasm_bindgen(typescript_type = "MoleculePropertyError")]
    pub type SdfReadErrorMoleculeProperty;
}
fn sdf_read_error_cause(e: &ck::SdfReadError) -> Option<SdfReadErrorCause> {
    match e {
        ck::SdfReadError::StereoGroup(v) => {
            Some(crate::model_errors::StereoGroupError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfReadError::MolPost(v) => {
            Some(crate::mol_post_errors::MolPostError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfReadError::Topology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SdfReadError::Coordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SdfReadError::QueryGraph(v) => {
            Some(crate::query_model_errors::QueryGraphError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfReadError::AtomProperty(v) => {
            Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfReadError::BondValue(v) => {
            Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfReadError::MoleculeProperty(v) => {
            Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfReadError::Record { source, .. } => {
            Some(SdfReadError::from_rust(source.as_ref().clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(SdfReadError,SdfReadErrorKind,ck::SdfReadError,{
StereoGroup=>ck::SdfReadError::StereoGroup(_),
MolPost=>ck::SdfReadError::MolPost(_),
Topology=>ck::SdfReadError::Topology(_),
Coordinates=>ck::SdfReadError::Coordinates(_),
QueryGraph=>ck::SdfReadError::QueryGraph(_),
AtomProperty=>ck::SdfReadError::AtomProperty(_),
BondValue=>ck::SdfReadError::BondValue(_),
MoleculeProperty=>ck::SdfReadError::MoleculeProperty(_),
QueryRecord=>ck::SdfReadError::QueryRecord,
Empty=>ck::SdfReadError::Empty,
Counts=>ck::SdfReadError::Counts,
MissingCapability=>ck::SdfReadError::MissingCapability(_),
Unsupported=>ck::SdfReadError::Unsupported(_),
Parse=>ck::SdfReadError::Parse(_),
Field=>ck::SdfReadError::Field{..},
PropertyListCount=>ck::SdfReadError::PropertyListCount{..},
RecordIndexOutOfRange=>ck::SdfReadError::RecordIndexOutOfRange{..},
Record=>ck::SdfReadError::Record{..},
},{
#[wasm_bindgen(getter,js_name=stereoGroup)] pub fn stereo_group(&self)->Option<SdfReadErrorStereoGroup>{match self.inner.as_ref(){ck::SdfReadError::StereoGroup(v)=>Some(crate::model_errors::StereoGroupError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=molPost)] pub fn mol_post(&self)->Option<SdfReadErrorMolPost>{match self.inner.as_ref(){ck::SdfReadError::MolPost(v)=>Some(crate::mol_post_errors::MolPostError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=topology)] pub fn topology(&self)->Option<SdfReadErrorTopology>{match self.inner.as_ref(){ck::SdfReadError::Topology(v)=>Some(crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=coordinates)] pub fn coordinates(&self)->Option<SdfReadErrorCoordinates>{match self.inner.as_ref(){ck::SdfReadError::Coordinates(v)=>Some(crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=queryGraph)] pub fn query_graph(&self)->Option<SdfReadErrorQueryGraph>{match self.inner.as_ref(){ck::SdfReadError::QueryGraph(v)=>Some(crate::query_model_errors::QueryGraphError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atomProperty)] pub fn atom_property(&self)->Option<SdfReadErrorAtomProperty>{match self.inner.as_ref(){ck::SdfReadError::AtomProperty(v)=>Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=bondValue)] pub fn bond_value(&self)->Option<SdfReadErrorBondValue>{match self.inner.as_ref(){ck::SdfReadError::BondValue(v)=>Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=moleculeProperty)] pub fn molecule_property(&self)->Option<SdfReadErrorMoleculeProperty>{match self.inner.as_ref(){ck::SdfReadError::MoleculeProperty(v)=>Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=missingCapability)] pub fn missing_capability(&self)->Option<String>{match self.inner.as_ref(){ck::SdfReadError::MissingCapability(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=unsupported)] pub fn unsupported(&self)->Option<String>{match self.inner.as_ref(){ck::SdfReadError::Unsupported(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=parse)] pub fn parse(&self)->Option<String>{match self.inner.as_ref(){ck::SdfReadError::Parse(v)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=kind)] pub fn kind(&self)->Option<String>{match self.inner.as_ref(){
ck::SdfReadError::Field{kind,..}=>Some(kind.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=line)] pub fn line(&self)->Option<usize>{match self.inner.as_ref(){
ck::SdfReadError::Field{line,..}=>Some(*line),
_=>None}}
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::SdfReadError::Field{value,..}=>Some(value.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=target)] pub fn target(&self)->Option<String>{match self.inner.as_ref(){
ck::SdfReadError::PropertyListCount{target,..}=>Some(target.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=name)] pub fn name(&self)->Option<String>{match self.inner.as_ref(){
ck::SdfReadError::PropertyListCount{name,..}=>Some(name.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=actual)] pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){
ck::SdfReadError::PropertyListCount{actual,..}=>Some(*actual),
_=>None}}
#[wasm_bindgen(getter,js_name=expected)] pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){
ck::SdfReadError::PropertyListCount{expected,..}=>Some(*expected),
_=>None}}
#[wasm_bindgen(getter,js_name=index)] pub fn index(&self)->Option<usize>{match self.inner.as_ref(){
ck::SdfReadError::RecordIndexOutOfRange{index,..}=>Some(*index),
ck::SdfReadError::Record{index,..}=>Some(*index),
_=>None}}
#[wasm_bindgen(getter,js_name=recordCount)] pub fn record_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::SdfReadError::RecordIndexOutOfRange{record_count,..}=>Some(*record_count),
_=>None}}
#[wasm_bindgen(getter,js_name=byteOffset)] pub fn byte_offset(&self)->Option<u64>{match self.inner.as_ref(){
ck::SdfReadError::Record{byte_offset,..}=>Some(*byte_offset),
_=>None}}
#[wasm_bindgen(getter,js_name=lineOffset)] pub fn line_offset(&self)->Option<usize>{match self.inner.as_ref(){
ck::SdfReadError::Record{line_offset,..}=>Some(*line_offset),
_=>None}}
#[wasm_bindgen(getter)] pub fn source(&self)->Option<SdfReadErrorCause>{match self.inner.as_ref(){ck::SdfReadError::Record{source,..}=>Some(SdfReadError::from_rust(source.as_ref().clone()).unchecked_into()),_=>None}}
},source=SdfReadErrorCause,cause=sdf_read_error_cause);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "SdfReadError | MolPostError | OperationError | QueryGraphError"
    )]
    #[derive(Clone)]
    pub type SdfErrorCause;
    #[wasm_bindgen(typescript_type = "SdfReadError")]
    pub type SdfErrorRead;
    #[wasm_bindgen(typescript_type = "MolPostError")]
    pub type SdfErrorPost;
    #[wasm_bindgen(typescript_type = "OperationError")]
    pub type SdfErrorConstruction;
    #[wasm_bindgen(typescript_type = "QueryGraphError")]
    pub type SdfErrorQueryGraph;
}
fn sdf_error_cause(e: &ck::SdfError) -> Option<SdfErrorCause> {
    match e {
        ck::SdfError::Read(v) => {
            Some(crate::sdf_errors::SdfReadError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfError::Post(v) => {
            Some(crate::mol_post_errors::MolPostError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfError::Construction(v) => {
            Some(crate::operation_errors::OperationError::from_rust(v.clone()).unchecked_into())
        }
        ck::SdfError::QueryGraph(v) => {
            Some(crate::query_model_errors::QueryGraphError::from_rust(v.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(SdfError,SdfErrorKind,ck::SdfError,{
Read=>ck::SdfError::Read(_),
Post=>ck::SdfError::Post(_),
Construction=>ck::SdfError::Construction(_),
QueryGraph=>ck::SdfError::QueryGraph(_),
QueryRecord=>ck::SdfError::QueryRecord,
WrongGraphKind=>ck::SdfError::WrongGraphKind{..},
},{
#[wasm_bindgen(getter,js_name=read)] pub fn read(&self)->Option<SdfErrorRead>{match self.inner.as_ref(){ck::SdfError::Read(v)=>Some(crate::sdf_errors::SdfReadError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=post)] pub fn post(&self)->Option<SdfErrorPost>{match self.inner.as_ref(){ck::SdfError::Post(v)=>Some(crate::mol_post_errors::MolPostError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=construction)] pub fn construction(&self)->Option<SdfErrorConstruction>{match self.inner.as_ref(){ck::SdfError::Construction(v)=>Some(crate::operation_errors::OperationError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=queryGraph)] pub fn query_graph(&self)->Option<SdfErrorQueryGraph>{match self.inner.as_ref(){ck::SdfError::QueryGraph(v)=>Some(crate::query_model_errors::QueryGraphError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=expected)] pub fn expected(&self)->Option<String>{match self.inner.as_ref(){
ck::SdfError::WrongGraphKind{expected,..}=>Some(expected.to_string()),
_=>None}}
#[wasm_bindgen(getter,js_name=actual)] pub fn actual(&self)->Option<String>{match self.inner.as_ref(){
ck::SdfError::WrongGraphKind{actual,..}=>Some(actual.to_string()),
_=>None}}
},source=SdfErrorCause,cause=sdf_error_cause);
