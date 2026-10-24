//! XYZ failures preserve numeric-parser and model sources.
use crate::errors::{enum_error, error_projection, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

error_projection!(ParseFloatError, std::num::ParseFloatError, {});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "ParseFloatError | TopologyValidationError | CoordinateValidationError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type XyzReadCause;
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | CoordinateValidationError | PropertyStringError"
    )]
    #[derive(Clone)]
    pub type XyzWriteCause;
}
fn read_cause(error: &ck::XyzReadError) -> Option<XyzReadCause> {
    match error {
        ck::XyzReadError::Coordinate { source, .. } => {
            Some(ParseFloatError::from_rust(source.clone()).unchecked_into())
        }
        ck::XyzReadError::Topology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::XyzReadError::Coordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::XyzReadError::MoleculeProperty(v) => {
            Some(crate::model_errors::MoleculePropertyError::from_rust(v.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(XyzReadError, XyzReadErrorKind, ck::XyzReadError, {
    EmptyBlock=>ck::XyzReadError::EmptyBlock, AtomCount=>ck::XyzReadError::AtomCount{..},
    UnexpectedEof=>ck::XyzReadError::UnexpectedEof, MissingCoordinates=>ck::XyzReadError::MissingCoordinates{..},
    Coordinate=>ck::XyzReadError::Coordinate{..}, AtomSymbol=>ck::XyzReadError::AtomSymbol{..},
    Topology=>ck::XyzReadError::Topology(_), Coordinates=>ck::XyzReadError::Coordinates(_), MoleculeProperty=>ck::XyzReadError::MoleculeProperty(_),
}, {
    #[wasm_bindgen(getter)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){ck::XyzReadError::AtomCount{value}|ck::XyzReadError::Coordinate{value,..}=>Some(value.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn line(&self)->Option<usize>{match self.inner.as_ref(){ck::XyzReadError::MissingCoordinates{line}|ck::XyzReadError::Coordinate{line,..}=>Some(*line),_=>None}}
    #[wasm_bindgen(getter)] pub fn message(&self)->Option<String>{match self.inner.as_ref(){ck::XyzReadError::AtomSymbol{message}=>Some(message.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<XyzReadCause>{read_cause(self.inner.as_ref())}
}, source=XyzReadCause,cause=read_cause);
fn write_cause(error: &ck::XyzWriteError) -> Option<XyzWriteCause> {
    match error {
        ck::XyzWriteError::Topology(v) => Some(
            crate::model_errors::TopologyValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::XyzWriteError::Coordinates(v) => Some(
            crate::model_errors::CoordinateValidationError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::XyzWriteError::Property(v) => {
            Some(crate::property_errors::PropertyStringError::from_rust(*v).unchecked_into())
        }
        ck::XyzWriteError::ConformerNotFound { .. } => None,
    }
}
enum_error!(XyzWriteError,XyzWriteErrorKind,ck::XyzWriteError,{
    Topology=>ck::XyzWriteError::Topology(_),Coordinates=>ck::XyzWriteError::Coordinates(_),Property=>ck::XyzWriteError::Property(_),ConformerNotFound=>ck::XyzWriteError::ConformerNotFound{..},
},{
    #[wasm_bindgen(getter)] pub fn id(&self)->Option<usize>{match self.inner.as_ref(){ck::XyzWriteError::ConformerNotFound{id}=>Some(*id),_=>None}}
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<XyzWriteCause>{write_cause(self.inner.as_ref())}
},source=XyzWriteCause,cause=write_cause);
