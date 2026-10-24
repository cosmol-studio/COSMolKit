//! Alignment errors retain source-transparent matches and thread ownership.
use crate::errors::{enum_error, error_projection, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ThreadCountError")]
    #[derive(Clone)]
    pub type ThreadSelectionCause;
    #[wasm_bindgen(
        typescript_type = "NonNullable<SubstructMatchError['cause'] | QueryGraphError['cause'] | SmartsParseError['cause']> | ThreadCountError"
    )]
    #[derive(Clone)]
    pub type AlignmentCause;
    #[wasm_bindgen(
        typescript_type = "CoordinateValidationError | SubstructMatchError | QueryGraphError | SmartsParseError | ThreadSelectionError"
    )]
    pub type AlignmentPayload;
}
fn thread_source(owner: &ReadonlyValue<ck::ThreadSelectionError>) -> Option<ThreadSelectionCause> {
    Some(
        crate::worker_errors::ThreadCountError::from_owner(owner.child(|e| {
            e.source()
                .expect("native thread selection source")
                .downcast_ref::<ck::ThreadCountError>()
                .expect("native thread selection source type")
        }))
        .unchecked_into(),
    )
}
error_projection!(
    ThreadSelectionError,
    ck::ThreadSelectionError,
    {
        #[wasm_bindgen(getter)]
        pub fn source(&self) -> Option<ThreadSelectionCause> {
            thread_source(&self.inner)
        }
    },
    source = ThreadSelectionCause,
    cause = thread_source
);
fn source(owner: &ReadonlyValue<ck::AlignmentError>) -> Option<AlignmentCause> {
    use ck::AlignmentError as E;
    match owner.as_ref() {
        E::Matching(v) => crate::match_errors::match_source(v).map(JsCast::unchecked_into),
        E::QueryGraph(v) => crate::query_model_errors::cause(v).map(JsCast::unchecked_into),
        E::QueryParse(v) => crate::smarts_parse_errors::parse_source(v).map(JsCast::unchecked_into),
        E::ThreadSelection(_) => thread_source(&owner.child(|e| match e {
            E::ThreadSelection(v) => v,
            _ => unreachable!("immutable selected thread variant"),
        }))
        .map(JsCast::unchecked_into),
        _ => None,
    }
}
enum_error!(AlignmentError,AlignmentErrorKind,ck::AlignmentError,{
    InvalidCoordinates=>ck::AlignmentError::InvalidCoordinates(_),Matching=>ck::AlignmentError::Matching(_),QueryGraph=>ck::AlignmentError::QueryGraph(_),QueryParse=>ck::AlignmentError::QueryParse(_),ThreadSelection=>ck::AlignmentError::ThreadSelection(_),NoConformers=>ck::AlignmentError::NoConformers,ConformerNotFound=>ck::AlignmentError::ConformerNotFound{..},EmptyAtomMap=>ck::AlignmentError::EmptyAtomMap,ProbeAtomOutOfRange=>ck::AlignmentError::ProbeAtomOutOfRange{..},ReferenceAtomOutOfRange=>ck::AlignmentError::ReferenceAtomOutOfRange{..},WeightCountMismatch=>ck::AlignmentError::WeightCountMismatch{..},NonPositiveWeight=>ck::AlignmentError::NonPositiveWeight{..},NoSubstructureMatch=>ck::AlignmentError::NoSubstructureMatch,TerminalGroupSymmetrization=>ck::AlignmentError::TerminalGroupSymmetrization{..},NumericalPrecondition=>ck::AlignmentError::NumericalPrecondition{..},WorkerTerminated=>ck::AlignmentError::WorkerTerminated,
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<AlignmentPayload>{use ck::AlignmentError as E;Some(match self.inner.as_ref(){E::InvalidCoordinates(v)=>crate::model_errors::CoordinateValidationError::from_rust(v.clone()),E::Matching(v)=>crate::match_errors::SubstructMatchError::from_rust(v.clone()),E::QueryGraph(v)=>crate::query_model_errors::QueryGraphError::from_rust(v.clone()),E::QueryParse(v)=>crate::smarts_parse_errors::SmartsParseError::from_rust(v.clone()),E::ThreadSelection(v)=>ThreadSelectionError::from_rust(v.clone()),_=>return None}.unchecked_into())}
    #[wasm_bindgen(getter)] pub fn id(&self)->Option<i32>{match self.inner.as_ref(){ck::AlignmentError::ConformerNotFound{id}=>Some(*id),_=>None}}
    #[wasm_bindgen(getter)] pub fn index(&self)->Option<usize>{use ck::AlignmentError as E;match self.inner.as_ref(){E::ProbeAtomOutOfRange{index,..}|E::ReferenceAtomOutOfRange{index,..}|E::NonPositiveWeight{index}=>Some(*index),_=>None}}
    #[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{use ck::AlignmentError as E;match self.inner.as_ref(){E::ProbeAtomOutOfRange{atom_count,..}|E::ReferenceAtomOutOfRange{atom_count,..}=>Some(*atom_count),_=>None}}
    #[wasm_bindgen(getter,js_name=mapLen)] pub fn map_len(&self)->Option<usize>{match self.inner.as_ref(){ck::AlignmentError::WeightCountMismatch{map_len,..}=>Some(*map_len),_=>None}}
    #[wasm_bindgen(getter,js_name=weightLen)] pub fn weight_len(&self)->Option<usize>{match self.inner.as_ref(){ck::AlignmentError::WeightCountMismatch{weight_len,..}=>Some(*weight_len),_=>None}}
    #[wasm_bindgen(getter)] pub fn message(&self)->Option<String>{use ck::AlignmentError as E;match self.inner.as_ref(){E::TerminalGroupSymmetrization{message}|E::NumericalPrecondition{message}=>Some((*message).into()),_=>None}}
},source=AlignmentCause,cause=source);
