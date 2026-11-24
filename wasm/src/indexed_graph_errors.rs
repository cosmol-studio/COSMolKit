use crate::errors::{enum_error, exception};
use ck::IndexedGraphError as E;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyEditError | TopologyValidationError | AtomPathSubsetCopyError | PropertyUIntReadError | PropertyStringError | ValenceError"
    )]
    #[derive(Clone)]
    pub type GraphCause;
}
fn cause(e: &E) -> Option<GraphCause> {
    Some(
        match e {
            E::TopologyEdit(v) => crate::topology_errors::TopologyEditError::from_rust(v.clone()),
            E::Topology(v) => crate::model_errors::TopologyValidationError::from_rust(v.clone()),
            E::Subset(v) => crate::fragment_errors::AtomPathSubsetCopyError::from_rust(v.clone()),
            E::NumericProperty(v) => {
                crate::property_errors::PropertyUIntReadError::from_rust(v.clone())
            }
            E::TextProperty(v) => crate::property_errors::PropertyStringError::from_rust(v.clone()),
            E::Valence(v) => crate::valence_errors::ValenceError::from_rust(v.clone()),
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(IndexedGraphError,IndexedGraphErrorKind,E,{
    EmptyAtomOrder=>E::EmptyAtomOrder,InvalidIndexProperty=>E::InvalidIndexProperty{..},InvalidStereoNeighbor=>E::InvalidStereoNeighbor{..},CoordinateKey=>E::CoordinateKey{..},
    TopologyEdit=>E::TopologyEdit(_),Topology=>E::Topology(_),Subset=>E::Subset(_),NumericProperty=>E::NumericProperty(_),TextProperty=>E::TextProperty(_),Valence=>E::Valence(_),
},{
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{match *self.inner {E::InvalidIndexProperty{atom,..}|E::InvalidStereoNeighbor{atom,..}=>Some(atom.index()),_=>None}}
    #[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{match *self.inner {E::InvalidStereoNeighbor{bond,..}=>Some(bond.index()),_=>None}}
    #[wasm_bindgen(getter)] pub fn property(&self)->Option<String>{match *self.inner {E::InvalidIndexProperty{property,..}=>Some(property.into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn dimension(&self)->Option<u8>{match *self.inner {E::CoordinateKey{dimension,..}=>Some(dimension),_=>None}}
    #[wasm_bindgen(getter)] pub fn id(&self)->Option<usize>{match *self.inner {E::CoordinateKey{id,..}=>Some(id),_=>None}}
    #[wasm_bindgen(getter)] pub fn reason(&self)->Option<String>{match *self.inner {E::CoordinateKey{reason,..}=>Some(reason.into()),_=>None}}
},source=GraphCause,cause=cause);
