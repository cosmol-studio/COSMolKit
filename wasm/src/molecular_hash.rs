//! Native hash and caller-provided rank queries retain exact u64 results.
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fields::Field,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CipRankError | TopologyValidationError")]
    #[derive(Clone)]
    pub type HashCause;
}
fn cause(e: &ck::MoleculeHashError) -> Option<HashCause> {
    match e {
        ck::MoleculeHashError::CipRanks(v) => {
            Some(crate::cip_rank_errors::CipRankError::from_rust(v.clone()).unchecked_into())
        }
        ck::MoleculeHashError::InvalidTopology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(MoleculeHashError,MoleculeHashErrorKind,ck::MoleculeHashError,{
 EmptyMolecule=>ck::MoleculeHashError::EmptyMolecule,MissingPreparedValence=>ck::MoleculeHashError::MissingPreparedValence,CipRanks=>ck::MoleculeHashError::CipRanks(_),InvalidTopology=>ck::MoleculeHashError::InvalidTopology(_),RankCount=>ck::MoleculeHashError::RankCount{..}
},{
 #[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{match *self.inner{ck::MoleculeHashError::RankCount{actual,..}=>Some(actual),_=>None}}
 #[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match *self.inner{ck::MoleculeHashError::RankCount{atom_count,..}=>Some(atom_count),_=>None}}
 #[wasm_bindgen(getter)]pub fn source(&self)->Option<HashCause>{cause(&self.inner)}
},source=HashCause,cause=cause);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=molecularHash)]
    pub fn molecular_hash(&self) -> Result<u64, JsValue> {
        self.0
            .molecular_hash()
            .map_err(MoleculeHashError::from_rust)
    }
    #[wasm_bindgen(js_name=molecularHashWithRanks)]
    pub fn molecular_hash_with_ranks(
        &self,
        ranks: <Vec<u32> as Field>::Input,
    ) -> Result<u64, JsValue> {
        self.0
            .molecular_hash_with_ranks(&Vec::<u32>::input(ranks)?)
            .map_err(MoleculeHashError::from_rust)
    }
}
