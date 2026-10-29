//! MOL writer payloads and source-transparent dependency projections.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {}
fn source(e: &ck::MolWriteError) -> Option<MolWriteCause> {
    use ck::MolWriteError as E;
    match e {
        #[cfg(any(feature = "cap-depict", feature = "cap-fingerprints"))]
        E::Depict(v) => crate::depict_errors::cause(v).map(JsCast::unchecked_into),
        E::QueryGraph(v) => crate::query_model_errors::cause(v).map(JsCast::unchecked_into),
        E::QueryState(v) => crate::query_errors::cause(v).map(JsCast::unchecked_into),
        #[cfg(feature = "cap-search")]
        E::QuerySmarts(v) => crate::smarts_errors::write_source(v).map(JsCast::unchecked_into),
        E::Valence(v) => {
            Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into())
        }
        E::Kekulize(v) => crate::kekulize_errors::cause(v).map(JsCast::unchecked_into),
        E::Atropisomer(v) => crate::atropisomer_errors::cause(v).map(JsCast::unchecked_into),
        E::Wedge(v) => crate::wedge_errors::cause(v).map(JsCast::unchecked_into),
        E::Topology(v) => crate::model_errors::topology_cause(v).map(JsCast::unchecked_into),
        E::Detached(v) => crate::io_write_errors::sdf_source(v).map(JsCast::unchecked_into),
        E::QueryAtom(_)
        | E::Coordinates(_)
        | E::Property(_)
        | E::UnsignedProperty(_)
        | E::MissingCapability(_)
        | E::UnsupportedSubset(_)
        | E::Value(_)
        | E::AmbiguousCoordinates { .. }
        | E::MissingCoordinate { .. }
        | E::CoordinateDimensionMismatch => None,
    }
}
enum_error!(MolWriteError,MolWriteErrorKind,ck::MolWriteError,{
    #[cfg(any(feature = "cap-depict", feature = "cap-fingerprints"))]
Depict=>ck::MolWriteError::Depict(_),
    MissingCapability=>ck::MolWriteError::MissingCapability(_),UnsupportedSubset=>ck::MolWriteError::UnsupportedSubset(_),Value=>ck::MolWriteError::Value(_),
    AmbiguousCoordinates=>ck::MolWriteError::AmbiguousCoordinates{..},MissingCoordinate=>ck::MolWriteError::MissingCoordinate{..},CoordinateDimensionMismatch=>ck::MolWriteError::CoordinateDimensionMismatch,
    QueryGraph=>ck::MolWriteError::QueryGraph(_),QueryAtom=>ck::MolWriteError::QueryAtom(_),QueryState=>ck::MolWriteError::QueryState(_),#[cfg(feature = "cap-search")]
QuerySmarts=>ck::MolWriteError::QuerySmarts(_),
    Valence=>ck::MolWriteError::Valence(_),Kekulize=>ck::MolWriteError::Kekulize(_),Atropisomer=>ck::MolWriteError::Atropisomer(_),Wedge=>ck::MolWriteError::Wedge(_),
    Topology=>ck::MolWriteError::Topology(_),Coordinates=>ck::MolWriteError::Coordinates(_),Property=>ck::MolWriteError::Property(_),UnsignedProperty=>ck::MolWriteError::UnsignedProperty(_),Detached=>ck::MolWriteError::Detached(_),
},{
    #[wasm_bindgen(getter,js_name=missingCapability)] pub fn missing_capability(&self)->Option<String>{match self.inner.as_ref(){ck::MolWriteError::MissingCapability(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=unsupportedSubset)] pub fn unsupported_subset(&self)->Option<String>{match self.inner.as_ref(){ck::MolWriteError::UnsupportedSubset(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){ck::MolWriteError::Value(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter,js_name=twoD)] pub fn two_d(&self)->Option<usize>{match self.inner.as_ref(){ck::MolWriteError::AmbiguousCoordinates{two_d,..}=>Some(*two_d),_=>None}}
    #[wasm_bindgen(getter,js_name=threeD)] pub fn three_d(&self)->Option<usize>{match self.inner.as_ref(){ck::MolWriteError::AmbiguousCoordinates{three_d,..}=>Some(*three_d),_=>None}}
    #[wasm_bindgen(getter)] pub fn id(&self)->Option<usize>{match self.inner.as_ref(){ck::MolWriteError::MissingCoordinate{id,..}=>Some(*id),_=>None}}
    #[wasm_bindgen(getter)] pub fn dimension(&self)->Option<crate::enums::CoordinateDimension>{match self.inner.as_ref(){ck::MolWriteError::MissingCoordinate{dimension,..}=>Some((*dimension).into()),_=>None}}
    // The payload is distinct from cause for transparent native variants.
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<MolWritePayload>{use ck::MolWriteError as E;Some(match self.inner.as_ref(){
#[cfg(any(feature = "cap-depict", feature = "cap-fingerprints"))]
        E::Depict(v)=>crate::depict_errors::Coordinate2DError::from_rust(v.clone()),
        E::QueryGraph(v)=>crate::query_model_errors::QueryGraphError::from_rust(v.clone()),
        E::QueryAtom(v)=>crate::query_model_errors::QueryAtomConversionError::from_rust(v.clone()),
        E::QueryState(v)=>crate::query_errors::QueryStateError::from_rust(v.clone()),
#[cfg(feature = "cap-search")]
        E::QuerySmarts(v)=>crate::smarts_errors::SmartsWriteError::from_rust(v.clone()),
        E::Valence(v)=>crate::valence_errors::ValenceError::from_rust(v.clone()),
        E::Kekulize(v)=>crate::kekulize_errors::KekulizeError::from_rust(v.clone()),
        E::Atropisomer(v)=>crate::atropisomer_errors::AtropisomerError::from_rust(v.clone()),
        E::Wedge(v)=>crate::wedge_errors::WedgeError::from_rust(v.clone()),
        E::Topology(v)=>crate::model_errors::TopologyValidationError::from_rust(v.clone()),
        E::Coordinates(v)=>crate::model_errors::CoordinateValidationError::from_rust(v.clone()),
        E::Property(v)=>crate::model_errors::PropertyValueError::from_rust(v.clone()),
        E::UnsignedProperty(v)=>crate::property_errors::PropertyUIntReadError::from_rust(v.clone()),
        E::Detached(v)=>crate::io_write_errors::SdfWriteError::from_rust(v.clone()),
        _=>return None,
    }.unchecked_into())}
},source=MolWriteCause,cause=source);

cosmolkit_tsify_macros::type_union! {
    MolWriteCause;
    "NonNullable<QueryGraphError['cause'] | QueryStateError['cause'] | KekulizeError['cause'] | AtropisomerError['cause'] | WedgeError['cause'] | TopologyValidationError['cause'] | SdfWriteError['cause']> | ValenceError";
    #[cfg(feature = "cap-search")]
    "NonNullable<SmartsWriteError['cause']>";
    #[cfg(any(feature = "cap-depict", feature = "cap-fingerprints"))]
    "NonNullable<Coordinate2DError['cause']>";
}

cosmolkit_tsify_macros::type_union! {
    MolWritePayload;
    "QueryGraphError | QueryAtomConversionError | QueryStateError | ValenceError | KekulizeError | AtropisomerError | WedgeError | TopologyValidationError | CoordinateValidationError | PropertyValueError | PropertyUIntReadError | SdfWriteError";
    #[cfg(feature = "cap-search")]
    "SmartsWriteError";
    #[cfg(any(feature = "cap-depict", feature = "cap-fingerprints"))]
    "Coordinate2DError";
}
