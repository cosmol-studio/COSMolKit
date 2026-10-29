use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "NonNullable<SubstructMatchError['cause']>")]
    #[derive(Clone)]
    pub type MatchCause;
    #[wasm_bindgen(typescript_type = "SubstructMatchError")]
    pub type MatchPayload;
}
fn cause(error: &ck::MatchError) -> Option<MatchCause> {
    match error {
        ck::MatchError::Substruct(v) => {
            crate::match_errors::match_source(v).map(|v| v.unchecked_into())
        }
        _ => None,
    }
}
enum_error!(MatchError,MatchErrorKind,ck::MatchError,{
    InvalidQuery=>ck::MatchError::InvalidQuery(_),InvalidTarget=>ck::MatchError::InvalidTarget(_),UnsupportedAtomPredicate=>ck::MatchError::UnsupportedAtomPredicate(_),UnsupportedBondPredicate=>ck::MatchError::UnsupportedBondPredicate(_),Substruct=>ck::MatchError::Substruct(_)
},{
    #[wasm_bindgen(getter)]pub fn message(&self)->Option<String>{match self.inner.as_ref(){ck::MatchError::InvalidQuery(v)|ck::MatchError::InvalidTarget(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter)]pub fn predicate(&self)->Option<String>{match self.inner.as_ref(){ck::MatchError::UnsupportedAtomPredicate(v)|ck::MatchError::UnsupportedBondPredicate(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn source(&self)->Option<MatchPayload>{match self.inner.as_ref(){ck::MatchError::Substruct(v)=>Some(crate::match_errors::SubstructMatchError::from_rust(v.clone()).unchecked_into()),_=>None}}
},source=MatchCause,cause=cause);
