//! Exact native torsion path scoring; bigint preserves the u64 result.
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fields::Field,
    input,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "TopologyValidationError | AtomCodeError | FingerprintError")]
    #[derive(Clone)]
    pub type PathScoreCause;
}
fn cause(e: &ck::TopologicalTorsionPathScoreError) -> Option<PathScoreCause> {
    use ck::TopologicalTorsionPathScoreError as E;
    Some(
        match e {
            E::InvalidTopology(v) => {
                crate::model_errors::TopologyValidationError::from_rust(v.clone())
            }
            E::AtomCode(v) => crate::fingerprint_errors::AtomCodeError::from_rust(v.clone()),
            E::PackedCode(v) => crate::fingerprint_errors::FingerprintError::from_rust(v.clone()),
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(TopologicalTorsionPathScoreError,TopologicalTorsionPathScoreErrorKind,ck::TopologicalTorsionPathScoreError,{
ZeroSize=>ck::TopologicalTorsionPathScoreError::ZeroSize,ShortPath=>ck::TopologicalTorsionPathScoreError::ShortPath{..},ShortAtomCodes=>ck::TopologicalTorsionPathScoreError::ShortAtomCodes{..},AtomIndexOutOfRange=>ck::TopologicalTorsionPathScoreError::AtomIndexOutOfRange{..},AtomCodeUnderflow=>ck::TopologicalTorsionPathScoreError::AtomCodeUnderflow{..},InvalidTopology=>ck::TopologicalTorsionPathScoreError::InvalidTopology(_),AtomCode=>ck::TopologicalTorsionPathScoreError::AtomCode(_),PackedCode=>ck::TopologicalTorsionPathScoreError::PackedCode(_)
},{
#[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{use ck::TopologicalTorsionPathScoreError as E;match *self.inner{E::ShortPath{actual,..}|E::ShortAtomCodes{actual,..}=>Some(actual),_=>None}}
#[wasm_bindgen(getter)]pub fn required(&self)->Option<usize>{use ck::TopologicalTorsionPathScoreError as E;match *self.inner{E::ShortPath{required,..}|E::ShortAtomCodes{required,..}=>Some(required),_=>None}}
#[wasm_bindgen(getter)]pub fn index(&self)->Option<usize>{use ck::TopologicalTorsionPathScoreError as E;match *self.inner{E::AtomIndexOutOfRange{index,..}|E::AtomCodeUnderflow{index,..}=>Some(index),_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match *self.inner{ck::TopologicalTorsionPathScoreError::AtomIndexOutOfRange{atom_count,..}=>Some(atom_count),_=>None}}
#[wasm_bindgen(getter)]pub fn code(&self)->Option<u32>{match *self.inner{ck::TopologicalTorsionPathScoreError::AtomCodeUnderflow{code,..}=>Some(code),_=>None}}
#[wasm_bindgen(getter)]pub fn subtract(&self)->Option<u32>{match *self.inner{ck::TopologicalTorsionPathScoreError::AtomCodeUnderflow{subtract,..}=>Some(subtract),_=>None}}
#[wasm_bindgen(getter)]pub fn source(&self)->Option<PathScoreCause>{cause(&self.inner)}
},source=PathScoreCause,cause=cause);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=topologicalTorsionPathScore)]
    pub fn topological_torsion_path_score(
        &self,
        path: <Vec<usize> as Field>::Input,
        size: tsify::Ts<input::IndexInput>,
        atom_codes: Option<<Vec<u32> as Field>::Input>,
    ) -> Result<u64, JsValue> {
        self.0
            .topological_torsion_path_score(
                &Vec::<usize>::input(path)?,
                input::read(size)?,
                atom_codes.map(Vec::<u32>::input).transpose()?.as_deref(),
            )
            .map_err(TopologicalTorsionPathScoreError::from_rust)
    }
}
