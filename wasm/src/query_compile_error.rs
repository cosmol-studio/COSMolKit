use crate::errors::enum_error;
use cosmolkit as ck;
use wasm_bindgen::prelude::*;

enum_error!(QueryCompileError,QueryCompileErrorKind,ck::QueryCompileError,{
    InvalidGraph=>ck::QueryCompileError::InvalidGraph(_),
},{
    #[wasm_bindgen(getter,js_name=invalidGraph)] pub fn invalid_graph(&self)->String{match self.inner.as_ref(){ck::QueryCompileError::InvalidGraph(v)=>v.clone()}}
});
