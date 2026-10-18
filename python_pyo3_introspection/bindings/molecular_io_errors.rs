//! Full molecular IO error graph shared by scalar, file and batch projections.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::bio_mmcif_reader::BioMmcifReadError;
use crate::bio_molecule_conversion::BioMoleculeError;
use crate::bio_pdb_reader::BioPdbReadError;
use crate::canonical_errors::PropertyStringError;
use crate::chemical_values::BondOrder;
use crate::coordinate_2d_errors::{Coordinate2DError, CoordinateSource};
use crate::embedding_values::CoordinateDimension;
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::IoError;
use crate::io_write_errors::{PdbWriteError, XyzWriteError};
use crate::kekulize_errors::{KekulizeError, KekulizeSource};
use crate::model_edit_errors::{MappingValidationError, QueryAtomConversionError, QueryStateError};
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use crate::mol2_post_errors::Mol2PostError;
use crate::mol2_read_errors::Mol2ReadError;
use crate::operation_errors::OperationError;
use crate::property_numeric_errors::{
    PropertyIntReadError, PropertyUIntReadError, UIntLexicalReadError,
};
use crate::property_values::PropertyValueError;
use crate::query_model_errors::QueryGraphError;
use crate::retained_error::RetainedError;
use crate::sdf_records::SdfError;
use crate::smarts_writer::{SmartsWriteError, WriterSource};
use crate::smarts_writer_support_errors::{WedgeError, WedgeSource};
use crate::stereo_errors::{AtropSource, AtropisomerError};
use crate::string_errors::AvalonEngineError;
use crate::utf8_errors::FromUtf8Error;
use crate::valence_errors::ValenceError;
use crate::xyz_errors::XyzReadError;
use cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
use std::error::Error;

pub(crate) type SdfWriteSource = Source<
    Py<PropertyStringError>,
    Source<
        Py<PropertyUIntReadError>,
        Source<
            Py<TopologyValidationError>,
            Source<Py<CoordinateValidationError>, Py<UIntLexicalReadError>>,
        >,
    >,
>;
error_projection!(SdfWriteError,SdfWriteErrorKind,ck::SdfWriteError,"io",{
Property(_),UnsignedProperty(_),UnsignedPropertyOverflow{atom,bond,property:String,value:u32},CountLimit,
Topology(_),Coordinates(_),IntegerProperty(_),PropertyValue(_),BondOrder(_),Atom(_),
SGroupBondMembership{group,bond},SubstanceGroup(_),
},cause=|inner:&ck::SdfWriteError|Python::attach(|py|match inner{
ck::SdfWriteError::Property(v)=>Some(PropertyStringError::from_rust(*v)),
ck::SdfWriteError::UnsignedProperty(v)=>Some(PropertyUIntReadError::from_rust(v.clone())),
ck::SdfWriteError::Topology(v)=>Some(TopologyValidationError::from_rust(v.clone())),
ck::SdfWriteError::Coordinates(v)=>Some(CoordinateValidationError::from_rust(v.clone())),
ck::SdfWriteError::IntegerProperty(v)=>PropertyIntReadError::from_rust(v.clone()).cause(py),
ck::SdfWriteError::PropertyValue(_)|ck::SdfWriteError::UnsignedPropertyOverflow{..}|ck::SdfWriteError::CountLimit|ck::SdfWriteError::BondOrder(_)|ck::SdfWriteError::Atom(_)|ck::SdfWriteError::SGroupBondMembership{..}|ck::SdfWriteError::SubstanceGroup(_)=>None,
}),source=Option<SdfWriteSource>,methods={
#[getter]fn atom(&self)->Option<usize>{match self.inner{ck::SdfWriteError::UnsignedPropertyOverflow{atom,..}=>atom.map(ck::AtomId::index),_=>None}}
#[getter]fn bond(&self)->Option<usize>{match self.inner{ck::SdfWriteError::UnsignedPropertyOverflow{bond,..}=>bond.map(ck::BondId::index),ck::SdfWriteError::SGroupBondMembership{bond,..}=>Some(bond.index()),_=>None}}
#[getter]fn group(&self)->Option<usize>{match &self.inner{ck::SdfWriteError::SGroupBondMembership{group,..}=>Some(group.index()),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn integer_property(&self,py:Python<'_>)->PyResult<Option<Py<PropertyIntReadError>>>{match &self.inner{ck::SdfWriteError::IntegerProperty(v)=>Ok(Some(PropertyIntReadError::from_rust(v.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn property_value(&self,py:Python<'_>)->PyResult<Option<Py<PropertyValueError>>>{match &self.inner{ck::SdfWriteError::PropertyValue(v)=>Ok(Some(PropertyValueError::from_rust(v.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn bond_order(&self)->Option<BondOrder>{match self.inner{ck::SdfWriteError::BondOrder(v)=>Some(v.into()),_=>None}}
#[getter]fn atom_state(&self)->Option<&'static str>{match self.inner{ck::SdfWriteError::Atom(v)=>Some(v),_=>None}}
#[getter]fn substance_group(&self)->Option<&str>{match &self.inner{ck::SdfWriteError::SubstanceGroup(v)=>Some(v),_=>None}}
});

type TopologySource = Source<Py<StereoGroupError>, Py<TemplateAttachmentOrderError>>;
pub(crate) type MolWriteSource = Source<
    Py<ValenceError>,
    Source<
        WriterSource,
        Source<
            KekulizeSource,
            Source<
                AtropSource,
                Source<
                    WedgeSource,
                    Source<
                        CoordinateSource,
                        Source<
                            TopologySource,
                            Source<
                                Py<MappingValidationError>,
                                Source<Py<UIntLexicalReadError>, SdfWriteSource>,
                            >,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;
error_projection!(retained MolWriteError,MolWriteErrorKind,ck::MolWriteError,"io",{
MissingCapability(_),UnsupportedSubset(_),Value(_),AmbiguousCoordinates{two_d:usize,three_d:usize},
MissingCoordinate{dimension:CoordinateDimension=>|v:&ck::CoordinateDimension|(*v).into(),id:usize},
CoordinateDimensionMismatch,QueryGraph(_),QueryAtom(_),QueryState(_),QuerySmarts(_),Valence(_),Kekulize(_),
Atropisomer(_),Wedge(_),Depict(_),Topology(_),Coordinates(_),Property(_),UnsignedProperty(_),Detached(_),
},cause=|inner:&RetainedError<ck::MolWriteError>|Python::attach(|py|match inner.as_ref(){
ck::MolWriteError::Valence(v)=>Some(ValenceError::from_rust(v.clone())),
ck::MolWriteError::QueryGraph(v)=>QueryGraphError::from_rust(v.clone()).cause(py),
ck::MolWriteError::QueryAtom(v)=>QueryAtomConversionError::from_rust(v.clone()).cause(py),
ck::MolWriteError::QueryState(v)=>QueryStateError::from_rust(v.clone()).cause(py),
ck::MolWriteError::QuerySmarts(_)=>SmartsWriteError::from_ref(inner.node().project(|e|match e.downcast_ref::<ck::MolWriteError>().unwrap(){ck::MolWriteError::QuerySmarts(v)=>v,_=>unreachable!()})).cause(py),
ck::MolWriteError::Kekulize(v)=>KekulizeError::from_rust(v.clone()).cause(py),
ck::MolWriteError::Atropisomer(v)=>AtropisomerError::from_rust(v.clone()).cause(py),
ck::MolWriteError::Wedge(v)=>WedgeError::from_rust(v.clone()).cause(py),
ck::MolWriteError::Depict(_)=>Coordinate2DError::from_ref(inner.node().project(|e|match e.downcast_ref::<ck::MolWriteError>().unwrap(){ck::MolWriteError::Depict(v)=>v,_=>unreachable!()})).cause(py),
ck::MolWriteError::Topology(v)=>TopologyValidationError::from_rust(v.clone()).cause(py),
ck::MolWriteError::Coordinates(v)=>CoordinateValidationError::from_rust(v.clone()).cause(py),
ck::MolWriteError::Property(v)=>PropertyValueError::from_rust(v.clone()).cause(py),
ck::MolWriteError::UnsignedProperty(v)=>PropertyUIntReadError::from_rust(v.clone()).cause(py),
ck::MolWriteError::Detached(v)=>SdfWriteError::from_rust(v.clone()).cause(py),
ck::MolWriteError::MissingCapability(_)|ck::MolWriteError::UnsupportedSubset(_)|ck::MolWriteError::Value(_)|ck::MolWriteError::AmbiguousCoordinates{..}|ck::MolWriteError::MissingCoordinate{..}|ck::MolWriteError::CoordinateDimensionMismatch=>None,
}),source=Option<MolWriteSource>,methods={
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn query_graph(&self,py:Python<'_>)->PyResult<Option<Py<QueryGraphError>>>{match self.inner.as_ref(){ck::MolWriteError::QueryGraph(value)=>Ok(Some(QueryGraphError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn query_atom(&self,py:Python<'_>)->PyResult<Option<Py<QueryAtomConversionError>>>{match self.inner.as_ref(){ck::MolWriteError::QueryAtom(value)=>Ok(Some(QueryAtomConversionError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn query_state(&self,py:Python<'_>)->PyResult<Option<Py<QueryStateError>>>{match self.inner.as_ref(){ck::MolWriteError::QueryState(value)=>Ok(Some(QueryStateError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn query_smarts(&self,py:Python<'_>)->PyResult<Option<Py<SmartsWriteError>>>{match self.inner.as_ref(){ck::MolWriteError::QuerySmarts(value)=>Ok(Some(SmartsWriteError::from_ref(self.inner.node().project(|e|match e.downcast_ref::<ck::MolWriteError>().unwrap(){ck::MolWriteError::QuerySmarts(value)=>value,_=>unreachable!()})).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn valence(&self,py:Python<'_>)->PyResult<Option<Py<ValenceError>>>{match self.inner.as_ref(){ck::MolWriteError::Valence(value)=>Ok(Some(ValenceError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn kekulize(&self,py:Python<'_>)->PyResult<Option<Py<KekulizeError>>>{match self.inner.as_ref(){ck::MolWriteError::Kekulize(value)=>Ok(Some(KekulizeError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn atropisomer(&self,py:Python<'_>)->PyResult<Option<Py<AtropisomerError>>>{match self.inner.as_ref(){ck::MolWriteError::Atropisomer(value)=>Ok(Some(AtropisomerError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn wedge(&self,py:Python<'_>)->PyResult<Option<Py<WedgeError>>>{match self.inner.as_ref(){ck::MolWriteError::Wedge(value)=>Ok(Some(WedgeError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn depict(&self,py:Python<'_>)->PyResult<Option<Py<Coordinate2DError>>>{match self.inner.as_ref(){ck::MolWriteError::Depict(value)=>Ok(Some(Coordinate2DError::from_ref(self.inner.node().project(|e|match e.downcast_ref::<ck::MolWriteError>().unwrap(){ck::MolWriteError::Depict(value)=>value,_=>unreachable!()})).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn topology(&self,py:Python<'_>)->PyResult<Option<Py<TopologyValidationError>>>{match self.inner.as_ref(){ck::MolWriteError::Topology(value)=>Ok(Some(TopologyValidationError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn coordinates(&self,py:Python<'_>)->PyResult<Option<Py<CoordinateValidationError>>>{match self.inner.as_ref(){ck::MolWriteError::Coordinates(value)=>Ok(Some(CoordinateValidationError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn property_error(&self,py:Python<'_>)->PyResult<Option<Py<PropertyValueError>>>{match self.inner.as_ref(){ck::MolWriteError::Property(value)=>Ok(Some(PropertyValueError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn unsigned_property(&self,py:Python<'_>)->PyResult<Option<Py<PropertyUIntReadError>>>{match self.inner.as_ref(){ck::MolWriteError::UnsignedProperty(value)=>Ok(Some(PropertyUIntReadError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn detached(&self,py:Python<'_>)->PyResult<Option<Py<SdfWriteError>>>{match self.inner.as_ref(){ck::MolWriteError::Detached(value)=>Ok(Some(SdfWriteError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]fn capability(&self)->Option<&'static str>{match self.inner.as_ref(){ck::MolWriteError::MissingCapability(v)=>Some(*v),_=>None}}
#[getter]fn unsupported_subset(&self)->Option<&'static str>{match self.inner.as_ref(){ck::MolWriteError::UnsupportedSubset(v)=>Some(*v),_=>None}}
#[getter]fn value(&self)->Option<&str>{match self.inner.as_ref(){ck::MolWriteError::Value(v)=>Some(v),_=>None}}
});

type MolecularIoSource = Source<
    Py<IoError>,
    Source<
        Py<XyzReadError>,
        Source<
            Py<XyzWriteError>,
            Source<
                Py<OperationError>,
                Source<
                    Py<Mol2ReadError>,
                    Source<
                        Py<Mol2PostError>,
                        Source<
                            Py<MolWriteError>,
                            Source<
                                Py<PdbWriteError>,
                                Source<
                                    Py<BioPdbReadError>,
                                    Source<
                                        Py<BioMmcifReadError>,
                                        Source<
                                            Py<BioMoleculeError>,
                                            Source<Py<SdfError>, Py<FromUtf8Error>>,
                                        >,
                                    >,
                                >,
                            >,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;
error_projection!(retained MolecularIoError,MolecularIoErrorKind,ck::MolecularIoError,"io",{
Io{path,source},XyzRead(_),XyzWrite(_),Construction(_),Mol2Read(_),Mol2Post(_),MolWrite(_),PdbWrite(_),
PdbRead(_),MmcifRead(_),BioMolecule(_),Sdf(_),NoRecord{format},
OutputUtf8{format,source},Parameter{name:&'static str,detail:&'static str},
},cause=|inner:&RetainedError<ck::MolecularIoError>|{
let node=inner.node().source()?;
Some(match inner.as_ref(){
ck::MolecularIoError::Io{..}=>IoError::from_ref(node),
ck::MolecularIoError::XyzRead(v)=>XyzReadError::from_rust(v.clone()),
ck::MolecularIoError::XyzWrite(v)=>XyzWriteError::from_rust(v.clone()),
ck::MolecularIoError::Construction(_)=>OperationError::from_ref(node),
ck::MolecularIoError::Mol2Read(v)=>Mol2ReadError::from_rust(v.clone()),
ck::MolecularIoError::Mol2Post(_)=>Mol2PostError::from_ref(node),
ck::MolecularIoError::MolWrite(_)=>MolWriteError::from_ref(node),
ck::MolecularIoError::PdbWrite(v)=>PdbWriteError::from_rust(v.clone()),
ck::MolecularIoError::PdbRead(_)=>BioPdbReadError::from_ref(node),
ck::MolecularIoError::MmcifRead(_)=>BioMmcifReadError::from_ref(node),
ck::MolecularIoError::BioMolecule(_)=>BioMoleculeError::from_ref(node),
ck::MolecularIoError::Sdf(_)=>SdfError::from_ref(node),
ck::MolecularIoError::OutputUtf8{source,..}=>FromUtf8Error::from_rust(source.clone()),
ck::MolecularIoError::NoRecord{..}|ck::MolecularIoError::Parameter{..}=>unreachable!("source-free IO error"),
})},source=Option<MolecularIoSource>,methods={
#[getter]fn path(&self)->Option<std::path::PathBuf>{match self.inner.as_ref(){ck::MolecularIoError::Io{path,..}=>Some(path.clone()),_=>None}}
#[getter]fn format(&self)->Option<&'static str>{match self.inner.as_ref(){ck::MolecularIoError::NoRecord{format}|ck::MolecularIoError::OutputUtf8{format,..}=>Some(*format),_=>None}}
});

error_projection!(retained AvalonFingerprintError,AvalonFingerprintErrorKind,ck::AvalonFingerprintError,"fingerprint",{
Input(_),Engine(_),
},cause=|inner:&RetainedError<ck::AvalonFingerprintError>|{
let node=inner.node().source()?;
Some(match inner.as_ref(){ck::AvalonFingerprintError::Input(_)=>MolecularIoError::from_ref(node),ck::AvalonFingerprintError::Engine(v)=>AvalonEngineError::from_rust(v.clone())})
},source=Option<Source<Py<MolecularIoError>,Py<AvalonEngineError>>>);
