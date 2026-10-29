//! Detached writer diagnostics; transparent sources follow native Error::source.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "PropertyStringError | PropertyUIntReadError | TopologyValidationError | CoordinateValidationError"
    )]
    #[derive(Clone)]
    pub type SdfWriteCause;
    #[wasm_bindgen(
        typescript_type = "PropertyStringError | PropertyUIntReadError | TopologyValidationError | CoordinateValidationError | PropertyIntReadError | PropertyValueError"
    )]
    pub type SdfWritePayload;
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | CoordinateValidationError | KekulizeError"
    )]
    #[derive(Clone)]
    pub type PdbWriteCause;
    #[wasm_bindgen(typescript_type = "PropertyStringError")]
    pub type PropertyStringPayload;
}
pub(crate) fn sdf_source(e: &ck::SdfWriteError) -> Option<SdfWriteCause> {
    use ck::SdfWriteError as E;
    Some(
        match e {
            E::Property(v) => crate::property_errors::PropertyStringError::from_rust(*v),
            E::UnsignedProperty(v) => {
                crate::property_errors::PropertyUIntReadError::from_rust(v.clone())
            }
            E::Topology(v) => crate::model_errors::TopologyValidationError::from_rust(v.clone()),
            E::Coordinates(v) => {
                crate::model_errors::CoordinateValidationError::from_rust(v.clone())
            }
            // Both transparent property leaves have no native source.
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(SdfWriteError,SdfWriteErrorKind,ck::SdfWriteError,{
    Property=>ck::SdfWriteError::Property(_),UnsignedProperty=>ck::SdfWriteError::UnsignedProperty(_),
    UnsignedPropertyOverflow=>ck::SdfWriteError::UnsignedPropertyOverflow{..},CountLimit=>ck::SdfWriteError::CountLimit,
    Topology=>ck::SdfWriteError::Topology(_),Coordinates=>ck::SdfWriteError::Coordinates(_),IntegerProperty=>ck::SdfWriteError::IntegerProperty(_),
    PropertyValue=>ck::SdfWriteError::PropertyValue(_),BondOrder=>ck::SdfWriteError::BondOrder(_),Atom=>ck::SdfWriteError::Atom(_),
    SGroupBondMembership=>ck::SdfWriteError::SGroupBondMembership{..},SubstanceGroup=>ck::SdfWriteError::SubstanceGroup(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<SdfWritePayload>{use ck::SdfWriteError as E;match self.inner.as_ref(){E::IntegerProperty(v)=>Some(crate::property_errors::PropertyIntReadError::from_rust(v.clone()).unchecked_into()),E::PropertyValue(v)=>Some(crate::model_errors::PropertyValueError::from_rust(v.clone()).unchecked_into()),v=>sdf_source(v).map(|v|v.unchecked_into())}}
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::SdfWriteError::UnsignedPropertyOverflow{atom,..}=>atom.map(|v|v.index()),_=>None}}
    #[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){ck::SdfWriteError::UnsignedPropertyOverflow{bond,..}=>bond.map(|v|v.index()),ck::SdfWriteError::SGroupBondMembership{bond,..}=>Some(bond.index()),_=>None}}
    #[wasm_bindgen(getter)] pub fn group(&self)->Option<usize>{match self.inner.as_ref(){ck::SdfWriteError::SGroupBondMembership{group,..}=>Some(group.index()),_=>None}}
    #[wasm_bindgen(getter)] pub fn property(&self)->Option<String>{match self.inner.as_ref(){ck::SdfWriteError::UnsignedPropertyOverflow{property,..}=>Some(property.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn value(&self)->Option<u32>{match self.inner.as_ref(){ck::SdfWriteError::UnsignedPropertyOverflow{value,..}=>Some(*value),_=>None}}
    #[wasm_bindgen(getter,js_name=bondOrder)] pub fn bond_order(&self)->Option<crate::enums::BondOrder>{match self.inner.as_ref(){ck::SdfWriteError::BondOrder(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=atomState)] pub fn atom_state(&self)->Option<String>{match self.inner.as_ref(){ck::SdfWriteError::Atom(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=substanceGroup)] pub fn substance_group(&self)->Option<String>{match self.inner.as_ref(){ck::SdfWriteError::SubstanceGroup(v)=>Some(v.clone()),_=>None}}
},source=SdfWriteCause,cause=sdf_source);

fn pdb_source(e: &ck::PdbWriteError) -> Option<PdbWriteCause> {
    use ck::PdbWriteError as E;
    Some(
        match e {
            E::Topology(v) => crate::model_errors::TopologyValidationError::from_rust(v.clone()),
            E::Coordinates(v) => {
                crate::model_errors::CoordinateValidationError::from_rust(v.clone())
            }
            E::Kekulize(v) => crate::kekulize_errors::KekulizeError::from_rust(v.clone()),
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(PdbWriteError,PdbWriteErrorKind,ck::PdbWriteError,{
    Property=>ck::PdbWriteError::Property(_),ConformerNotFound=>ck::PdbWriteError::ConformerNotFound{..},Unsupported=>ck::PdbWriteError::Unsupported{..},
    Topology=>ck::PdbWriteError::Topology(_),Coordinates=>ck::PdbWriteError::Coordinates(_),AmbiguousCoordinates=>ck::PdbWriteError::AmbiguousCoordinates{..},
    MissingCoordinate=>ck::PdbWriteError::MissingCoordinate{..},ConflictingCoordinateSelection=>ck::PdbWriteError::ConflictingCoordinateSelection,Kekulize=>ck::PdbWriteError::Kekulize(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<PdbWriteCause>{pdb_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn property(&self)->Option<PropertyStringPayload>{match self.inner.as_ref(){ck::PdbWriteError::Property(v)=>Some(crate::property_errors::PropertyStringError::from_rust(*v).unchecked_into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn id(&self)->Option<usize>{match self.inner.as_ref(){ck::PdbWriteError::ConformerNotFound{id}|ck::PdbWriteError::MissingCoordinate{id,..}=>Some(*id),_=>None}}
    #[wasm_bindgen(getter)] pub fn feature(&self)->Option<String>{match self.inner.as_ref(){ck::PdbWriteError::Unsupported{feature}=>Some((*feature).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=twoD)] pub fn two_d(&self)->Option<usize>{match self.inner.as_ref(){ck::PdbWriteError::AmbiguousCoordinates{two_d,..}=>Some(*two_d),_=>None}}
    #[wasm_bindgen(getter,js_name=threeD)] pub fn three_d(&self)->Option<usize>{match self.inner.as_ref(){ck::PdbWriteError::AmbiguousCoordinates{three_d,..}=>Some(*three_d),_=>None}}
    #[wasm_bindgen(getter)] pub fn dimension(&self)->Option<crate::enums::CoordinateDimension>{match self.inner.as_ref(){ck::PdbWriteError::MissingCoordinate{dimension,..}=>Some((*dimension).into()),_=>None}}
},source=PdbWriteCause,cause=pdb_source);
