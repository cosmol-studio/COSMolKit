//! Prepared-query and substructure matcher error projections.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | RingFindingError | NonNullable<ValenceError['cause']>"
    )]
    #[derive(Clone)]
    pub type QueryMatchContextCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError | RingFindingError | ValenceError")]
    pub type QueryMatchContextPayload;
    #[wasm_bindgen(
        typescript_type = "PropertyIntReadError | NonNullable<StereoOrderError['cause'] | QueryMatchContextError['cause']>"
    )]
    #[derive(Clone)]
    pub type SubstructMatchCause;
    #[wasm_bindgen(
        typescript_type = "PeriodicTableError | StereoOrderError | PropertyStringError | PropertyIntReadError | QueryMatchContextError"
    )]
    pub type SubstructMatchPayload;
}
pub(crate) fn context_source(e: &ck::QueryMatchContextError) -> Option<QueryMatchContextCause> {
    use ck::QueryMatchContextError as E;
    match e {
        E::InvalidTopology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        E::Rings(v) => {
            Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into())
        }
        E::Valence(v) => crate::valence_errors::cause(v).map(JsCast::unchecked_into),
        _ => None,
    }
}
enum_error!(QueryMatchContextError,QueryMatchContextErrorKind,ck::QueryMatchContextError,{
    InvalidTopology=>ck::QueryMatchContextError::InvalidTopology(_),ValenceRows=>ck::QueryMatchContextError::ValenceRows{..},Valence=>ck::QueryMatchContextError::Valence(_),ValencePrecondition=>ck::QueryMatchContextError::ValencePrecondition{..},RingPrecondition=>ck::QueryMatchContextError::RingPrecondition{..},UninitializedRings=>ck::QueryMatchContextError::UninitializedRings,RingMembershipRows=>ck::QueryMatchContextError::RingMembershipRows{..},Rings=>ck::QueryMatchContextError::Rings(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<QueryMatchContextPayload>{use ck::QueryMatchContextError as E;Some(match self.inner.as_ref(){E::InvalidTopology(v)=>crate::model_errors::TopologyValidationError::from_rust(v.clone()),E::Rings(v)=>crate::ring_errors::RingFindingError::from_rust(v.clone()),E::Valence(v)=>crate::valence_errors::ValenceError::from_rust(v.clone()),_=>return None}.unchecked_into())}
    #[wasm_bindgen(getter)] pub fn field(&self)->Option<String>{use ck::QueryMatchContextError as E;match self.inner.as_ref(){E::ValenceRows{field,..}|E::RingMembershipRows{field,..}|E::ValencePrecondition{field,..}=>Some((*field).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn expected(&self)->Option<usize>{use ck::QueryMatchContextError as E;match self.inner.as_ref(){E::ValenceRows{expected,..}|E::RingMembershipRows{expected,..}=>Some(*expected),_=>None}}
    #[wasm_bindgen(getter)] pub fn actual(&self)->Option<usize>{use ck::QueryMatchContextError as E;match self.inner.as_ref(){E::ValenceRows{actual,..}|E::RingMembershipRows{actual,..}=>Some(*actual),_=>None}}
    #[wasm_bindgen(getter)] pub fn getter(&self)->Option<String>{use ck::QueryMatchContextError as E;match self.inner.as_ref(){E::ValencePrecondition{getter,..}|E::RingPrecondition{getter,..}=>Some((*getter).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{use ck::QueryMatchContextError as E;match self.inner.as_ref(){E::ValencePrecondition{atom,..}=>Some(*atom),E::RingPrecondition{atom,..}=>*atom,_=>None}}
    #[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){ck::QueryMatchContextError::RingPrecondition{bond,..}=>*bond,_=>None}}
},source=QueryMatchContextCause,cause=context_source);
pub(crate) fn match_source(e: &ck::SubstructMatchError) -> Option<SubstructMatchCause> {
    use ck::SubstructMatchError as E;
    match e {
        E::StereoOrder(v) => crate::stereo_order_errors::cause(v).map(JsCast::unchecked_into),
        E::PropertyInteger { source, .. } => Some(
            crate::property_errors::PropertyIntReadError::from_rust(source.clone())
                .unchecked_into(),
        ),
        E::QueryContext(v) => context_source(v).map(JsCast::unchecked_into),
        _ => None,
    }
}
enum_error!(SubstructMatchError,SubstructMatchErrorKind,ck::SubstructMatchError,{
    Unsupported=>ck::SubstructMatchError::Unsupported{..},FinalCheckMappingLength=>ck::SubstructMatchError::FinalCheckMappingLength{..},FinalCheckMappingIndex=>ck::SubstructMatchError::FinalCheckMappingIndex{..},FinalCheckInvariant=>ck::SubstructMatchError::FinalCheckInvariant{..},FinalCheckBondEndpoint=>ck::SubstructMatchError::FinalCheckBondEndpoint{..},FinalCheckMissingBond=>ck::SubstructMatchError::FinalCheckMissingBond{..},PeriodicTable=>ck::SubstructMatchError::PeriodicTable(_),StereoOrder=>ck::SubstructMatchError::StereoOrder(_),PropertyString=>ck::SubstructMatchError::PropertyString(_),PropertyInteger=>ck::SubstructMatchError::PropertyInteger{..},QueryContext=>ck::SubstructMatchError::QueryContext(_),UninitializedRecursiveQuery=>ck::SubstructMatchError::UninitializedRecursiveQuery,
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<SubstructMatchPayload>{use ck::SubstructMatchError as E;Some(match self.inner.as_ref(){E::PeriodicTable(v)=>crate::periodic_errors::PeriodicTableError::from_rust(v.clone()),E::StereoOrder(v)=>crate::stereo_order_errors::StereoOrderError::from_rust(v.clone()),E::PropertyString(v)=>crate::property_errors::PropertyStringError::from_rust(*v),E::PropertyInteger{source,..}=>crate::property_errors::PropertyIntReadError::from_rust(source.clone()),E::QueryContext(v)=>QueryMatchContextError::from_rust(v.clone()),_=>return None}.unchecked_into())}
    #[wasm_bindgen(getter)] pub fn branch(&self)->Option<String>{match self.inner.as_ref(){ck::SubstructMatchError::Unsupported{branch,..}=>Some((*branch).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=rdkitFunction)] pub fn rdkit_function(&self)->Option<String>{match self.inner.as_ref(){ck::SubstructMatchError::Unsupported{rdkit_function,..}=>Some((*rdkit_function).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=queryMapping)] pub fn query_mapping(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckMappingLength{query_mapping,..}=>Some(*query_mapping),_=>None}}
    #[wasm_bindgen(getter,js_name=targetMapping)] pub fn target_mapping(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckMappingLength{target_mapping,..}=>Some(*target_mapping),_=>None}}
    #[wasm_bindgen(getter,js_name=queryAtoms)] pub fn query_atoms(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckMappingLength{query_atoms,..}=>Some(*query_atoms),_=>None}}
    #[wasm_bindgen(getter)] pub fn side(&self)->Option<String>{use ck::SubstructMatchError as E;match self.inner.as_ref(){E::FinalCheckMappingIndex{side,..}|E::FinalCheckBondEndpoint{side,..}=>Some((*side).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn position(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckMappingIndex{position,..}=>Some(*position),_=>None}}
    #[wasm_bindgen(getter)] pub fn index(&self)->Option<usize>{use ck::SubstructMatchError as E;match self.inner.as_ref(){E::FinalCheckMappingIndex{index,..}|E::FinalCheckBondEndpoint{index,..}=>Some(*index),_=>None}}
    #[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{use ck::SubstructMatchError as E;match self.inner.as_ref(){E::FinalCheckMappingIndex{atom_count,..}|E::FinalCheckBondEndpoint{atom_count,..}=>Some(*atom_count),_=>None}}
    #[wasm_bindgen(getter)] pub fn invariant(&self)->Option<String>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckInvariant{invariant,..}=>Some((*invariant).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=queryAtom)] pub fn query_atom(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckInvariant{query_atom,..}=>Some(*query_atom),_=>None}}
    #[wasm_bindgen(getter)] pub fn endpoint(&self)->Option<String>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckBondEndpoint{endpoint,..}=>Some((*endpoint).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=queryBond)] pub fn query_bond(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckMissingBond{query_bond,..}=>Some(*query_bond),_=>None}}
    #[wasm_bindgen(getter)] pub fn begin(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckMissingBond{begin,..}=>Some(*begin),_=>None}}
    #[wasm_bindgen(getter)] pub fn end(&self)->Option<usize>{match self.inner.as_ref(){ck::SubstructMatchError::FinalCheckMissingBond{end,..}=>Some(*end),_=>None}}
    #[wasm_bindgen(getter)] pub fn property(&self)->Option<String>{match self.inner.as_ref(){ck::SubstructMatchError::PropertyInteger{property,..}=>Some((*property).into()),_=>None}}
},source=SubstructMatchCause,cause=match_source);
