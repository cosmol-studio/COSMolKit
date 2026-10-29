//! Exact native tautomer execution errors.
use crate::errors::enum_error;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "SmilesParseError | ValenceError | RingFindingError | KekulizeError | SanitizeError | LegacyStereoError | QueryCompileError | SubstructMatchError | QueryMatchContextError | CoordinateValidationError | TopologyValidationError | AtomPropertyError | MoleculePropertyError | BondValueError | PropertyStringError"
    )]
    #[derive(Clone)]
    pub type TautomerRunCause;
}
fn cause(e: &ck::TautomerRunError) -> Option<TautomerRunCause> {
    match e {
        ck::TautomerRunError::Smiles(v) => Some(
            crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::TautomerRunError::Valence(v) => {
            Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::Rings(v) => {
            Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::Kekulize(v) => {
            Some(crate::kekulize_errors::KekulizeError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::Sanitize(v) => {
            Some(crate::sanitize_errors::SanitizeError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::LegacyStereo(v) => Some(
            crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::TautomerRunError::QueryCompile(v) => Some(
            crate::query_compile_error::QueryCompileError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::TautomerRunError::Match(v) => {
            Some(crate::match_errors::SubstructMatchError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::MatchContext(v) => {
            Some(crate::match_errors::QueryMatchContextError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::Coordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::TautomerRunError::Topology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::TautomerRunError::AtomProperty(v) => {
            Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::MoleculeProperty(v) => {
            Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::BondValue(v) => {
            Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into())
        }
        ck::TautomerRunError::PropertyString(v) => {
            Some(crate::property_errors::PropertyStringError::from_rust(v.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(TautomerRunError,TautomerRunErrorKind,ck::TautomerRunError,{
Smiles=>ck::TautomerRunError::Smiles(_),
Valence=>ck::TautomerRunError::Valence(_),
Rings=>ck::TautomerRunError::Rings(_),
Kekulize=>ck::TautomerRunError::Kekulize(_),
Sanitize=>ck::TautomerRunError::Sanitize(_),
LegacyStereo=>ck::TautomerRunError::LegacyStereo(_),
QueryCompile=>ck::TautomerRunError::QueryCompile(_),
Match=>ck::TautomerRunError::Match(_),
MatchContext=>ck::TautomerRunError::MatchContext(_),
Coordinates=>ck::TautomerRunError::Coordinates(_),
Topology=>ck::TautomerRunError::Topology(_),
AtomProperty=>ck::TautomerRunError::AtomProperty(_),
MoleculeProperty=>ck::TautomerRunError::MoleculeProperty(_),
BondValue=>ck::TautomerRunError::BondValue(_),
PropertyString=>ck::TautomerRunError::PropertyString(_),
AtomMappingCount=>ck::TautomerRunError::AtomMappingCount{..},
BondMappingCount=>ck::TautomerRunError::BondMappingCount{..},
AtomCountMismatch=>ck::TautomerRunError::AtomCountMismatch{..},
BondCountMismatch=>ck::TautomerRunError::BondCountMismatch{..},
AtomOutOfRange=>ck::TautomerRunError::AtomOutOfRange{..},
BondOutOfRange=>ck::TautomerRunError::BondOutOfRange{..},
EditCount=>ck::TautomerRunError::EditCount{..},
HydrogenCountOutOfRange=>ck::TautomerRunError::HydrogenCountOutOfRange{..},
FormalChargeOutOfRange=>ck::TautomerRunError::FormalChargeOutOfRange{..},
MissingMappedBond=>ck::TautomerRunError::MissingMappedBond{..},
MissingValence=>ck::TautomerRunError::MissingValence{..},
EmptyTransformMatch=>ck::TautomerRunError::EmptyTransformMatch,
NoCanonicalTautomer=>ck::TautomerRunError::NoCanonicalTautomer,
SourceCacheCommit=>ck::TautomerRunError::SourceCacheCommit,
ScoreCacheTargetMismatch=>ck::TautomerRunError::ScoreCacheTargetMismatch,
Control=>ck::TautomerRunError::Control(_),Callback=>ck::TautomerRunError::Callback(_)
},{
#[wasm_bindgen(getter,js_name=expected)] pub fn expected(&self)->Option<usize> { match self.inner.as_ref() {
ck::TautomerRunError::AtomMappingCount{expected,..}=>Some(*expected),
ck::TautomerRunError::BondMappingCount{expected,..}=>Some(*expected),
ck::TautomerRunError::AtomCountMismatch{expected,..}=>Some(*expected),
ck::TautomerRunError::BondCountMismatch{expected,..}=>Some(*expected),
ck::TautomerRunError::EditCount{expected,..}=>Some(*expected),
_=>None }}
#[wasm_bindgen(getter,js_name=actual)] pub fn actual(&self)->Option<usize> { match self.inner.as_ref() {
ck::TautomerRunError::AtomMappingCount{actual,..}=>Some(*actual),
ck::TautomerRunError::BondMappingCount{actual,..}=>Some(*actual),
ck::TautomerRunError::AtomCountMismatch{actual,..}=>Some(*actual),
ck::TautomerRunError::BondCountMismatch{actual,..}=>Some(*actual),
ck::TautomerRunError::EditCount{actual,..}=>Some(*actual),
_=>None }}
#[wasm_bindgen(getter,js_name=atom)] pub fn atom(&self)->Option<usize> { match self.inner.as_ref() {
ck::TautomerRunError::AtomOutOfRange{atom,..}=>Some(*atom),
ck::TautomerRunError::HydrogenCountOutOfRange{atom,..}=>Some(atom.index()),
ck::TautomerRunError::FormalChargeOutOfRange{atom,..}=>Some(atom.index()),
ck::TautomerRunError::MissingValence{atom,..}=>Some(atom.index()),
_=>None }}
#[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize> { match self.inner.as_ref() {
ck::TautomerRunError::AtomOutOfRange{atom_count,..}=>Some(*atom_count),
_=>None }}
#[wasm_bindgen(getter,js_name=bond)] pub fn bond(&self)->Option<usize> { match self.inner.as_ref() {
ck::TautomerRunError::BondOutOfRange{bond,..}=>Some(*bond),
ck::TautomerRunError::MissingMappedBond{bond,..}=>Some(*bond),
_=>None }}
#[wasm_bindgen(getter,js_name=bondCount)] pub fn bond_count(&self)->Option<usize> { match self.inner.as_ref() {
ck::TautomerRunError::BondOutOfRange{bond_count,..}=>Some(*bond_count),
_=>None }}
#[wasm_bindgen(getter,js_name=field)] pub fn field(&self)->Option<String> { match self.inner.as_ref() {
ck::TautomerRunError::EditCount{field,..}=>Some((*field).into()),
_=>None }}
#[wasm_bindgen(getter,js_name=count)] pub fn count(&self)->Option<u32> { match self.inner.as_ref() {
ck::TautomerRunError::HydrogenCountOutOfRange{count,..}=>Some(*count),
_=>None }}
#[wasm_bindgen(getter,js_name=charge)] pub fn charge(&self)->Option<i32> { match self.inner.as_ref() {
ck::TautomerRunError::FormalChargeOutOfRange{charge,..}=>Some(*charge),
_=>None }}
#[wasm_bindgen(getter)] pub fn control(&self)->Option<String>{match self.inner.as_ref(){ck::TautomerRunError::Control(v)=>Some(v.clone()),_=>None}}
#[wasm_bindgen(getter)] pub fn callback(&self)->Option<String>{match self.inner.as_ref(){ck::TautomerRunError::Callback(v)=>Some(v.clone()),_=>None}}
},source=TautomerRunCause,cause=cause);
