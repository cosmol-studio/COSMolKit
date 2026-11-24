use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{ErrorObject, Source, error_projection};
use crate::fragment_errors::AtomPathSubsetCopyError;
use crate::model_errors::TopologyValidationError;
use crate::property_numeric_errors::PropertyUIntReadError;
use crate::topology_edit_errors::TopologyEditError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;
type GraphCause = Source<
    ErrorObject<TopologyEditError>,
    Source<
        ErrorObject<TopologyValidationError>,
        Source<
            ErrorObject<AtomPathSubsetCopyError>,
            Source<
                ErrorObject<PropertyUIntReadError>,
                Source<ErrorObject<PropertyStringError>, ErrorObject<ValenceError>>,
            >,
        >,
    >,
>;

error_projection!(IndexedGraphError,IndexedGraphErrorKind,ck::IndexedGraphError,"indexed_graph",{
    EmptyAtomOrder, InvalidIndexProperty {atom, property:&'static str},
    InvalidStereoNeighbor {bond:usize => |v:&ck::BondId|v.index(),atom},
    CoordinateKey {dimension:u8,id:usize,reason:&'static str},
    TopologyEdit(_),Topology(_),Subset(_),NumericProperty(_),TextProperty(_),Valence(_),
},cause=|inner:&ck::IndexedGraphError| match inner {
    ck::IndexedGraphError::TopologyEdit(v)=>Some(TopologyEditError::from_rust(v.clone())),
    ck::IndexedGraphError::Topology(v)=>Some(TopologyValidationError::from_rust(v.clone())),
    ck::IndexedGraphError::Subset(v)=>Some(AtomPathSubsetCopyError::from_rust(v.clone())),
    ck::IndexedGraphError::NumericProperty(v)=>Some(PropertyUIntReadError::from_rust(v.clone())),
    ck::IndexedGraphError::TextProperty(v)=>Some(PropertyStringError::from_rust(v.clone())),
    ck::IndexedGraphError::Valence(v)=>Some(ValenceError::from_rust(v.clone())),
    _=>None,
},source=Option<GraphCause>,methods={
    #[getter]
    /// Atom carrying the invalid property or used as an invalid stereo neighbor.
    fn atom(&self)->Option<usize>{match self.inner {ck::IndexedGraphError::InvalidIndexProperty{atom,..}|ck::IndexedGraphError::InvalidStereoNeighbor{atom,..}=>Some(atom.index()),_=>None}}
});
