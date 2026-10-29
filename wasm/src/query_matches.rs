//! Compiled-query matches retain both native mappings without re-compilation.
use crate::{
    Molecule,
    fields::Field,
    match_error::MatchError,
    query_calls::{CompiledQuery, CompiledQueryInput},
};
use cosmolkit as ck;
use tsify::Tsify;
use wasm_bindgen::prelude::*;
#[derive(serde::Serialize, Tsify)]
#[serde(transparent)]
pub struct AtomPairs(Vec<(usize, usize)>);
#[wasm_bindgen]
pub struct MatchResult(pub(crate) ck::MatchResult);
crate::fields::value_field!(
    ck::MatchResult,
    MatchResult,
    MatchResultInput,
    "MatchResult"
);
crate::fields::value_array!(ck::MatchResult, MatchResults, "ReadonlyArray<MatchResult>");
#[wasm_bindgen]
impl MatchResult {
    #[wasm_bindgen(constructor)]
    pub fn new(
        atom_mapping: <Vec<usize> as Field>::Input,
        bond_mapping: <Vec<usize> as Field>::Input,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::MatchResult::new(
            Vec::<usize>::input(atom_mapping)?,
            Vec::<usize>::input(bond_mapping)?,
        )))
    }
    #[wasm_bindgen(getter,js_name=atomMapping)]
    pub fn atom_mapping(&self) -> Vec<usize> {
        self.0.atom_mapping.clone()
    }
    #[wasm_bindgen(getter,js_name=bondMapping)]
    pub fn bond_mapping(&self) -> Vec<usize> {
        self.0.bond_mapping.clone()
    }
    #[wasm_bindgen(js_name=atomPairs)]
    pub fn atom_pairs(&self) -> Result<tsify::Ts<AtomPairs>, JsValue> {
        AtomPairs(self.0.atom_pairs())
            .into_ts()
            .map_err(wasm_bindgen::JsError::from)
            .map_err(Into::into)
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=substructMatchesCompiled)]
    pub fn substruct_matches_compiled(
        &self,
        query: CompiledQueryInput,
    ) -> Result<MatchResults, JsValue> {
        CompiledQuery::with_ref(&query, |q| {
            Vec::<ck::MatchResult>::get(crate::storage::NativeField::owned(
                self.0
                    .substruct_matches_compiled(q)
                    .map_err(MatchError::from_rust)?,
            ))
        })
    }
}
