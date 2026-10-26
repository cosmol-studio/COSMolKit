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
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type SdfWriteSource = Source<
    crate::error_framework::ErrorObject<PropertyStringError>,
    Source<
        crate::error_framework::ErrorObject<PropertyUIntReadError>,
        Source<
            crate::error_framework::ErrorObject<TopologyValidationError>,
            Source<
                crate::error_framework::ErrorObject<CoordinateValidationError>,
                crate::error_framework::ErrorObject<UIntLexicalReadError>,
            >,
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
#[getter]#[doc = "Atom identifier implicated in the failure. Applies to SdfWriteError.UnsignedPropertyOverflow. None when this payload is absent."] fn atom(&self)->Option<usize>{match self.inner{ck::SdfWriteError::UnsignedPropertyOverflow{atom,..}=>atom.map(ck::AtomId::index),_=>None}}
#[getter]#[doc = "Bond identifier implicated in the failure. Applies to SdfWriteError.UnsignedPropertyOverflow, SdfWriteError.SGroupBondMembership. None when this payload is absent."] fn bond(&self)->Option<usize>{match self.inner{ck::SdfWriteError::UnsignedPropertyOverflow{bond,..}=>bond.map(ck::BondId::index),ck::SdfWriteError::SGroupBondMembership{bond,..}=>Some(bond.index()),_=>None}}
#[getter]#[doc = "Group identifier involved in the failure. Applies to SdfWriteError.SGroupBondMembership. None when this payload is absent."] fn group(&self)->Option<usize>{match &self.inner{ck::SdfWriteError::SGroupBondMembership{group,..}=>Some(group.index()),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured PropertyIntReadError payload of SdfWriteError.IntegerProperty; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn integer_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PropertyIntReadError>>>{match &self.inner{ck::SdfWriteError::IntegerProperty(v)=>Ok(Some(PropertyIntReadError::from_rust(v.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured PropertyValueError payload of SdfWriteError.PropertyValue; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn property_value(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PropertyValueError>>>{match &self.inner{ck::SdfWriteError::PropertyValue(v)=>Ok(Some(PropertyValueError::from_rust(v.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Bond order implicated in the failed chemistry check. Applies to SdfWriteError.BondOrder. None when this payload is absent."] fn bond_order(&self)->Option<BondOrder>{match self.inner{ck::SdfWriteError::BondOrder(v)=>Some(v.into()),_=>None}}
#[getter]#[doc="Atom state that prevented SDF writing; None for other variants."]fn atom_state(&self)->Option<&'static str>{match self.inner{ck::SdfWriteError::Atom(v)=>Some(v),_=>None}}
#[getter]#[doc="Substance-group diagnostic that prevented SDF writing; None for other variants."]fn substance_group(&self)->Option<&str>{match &self.inner{ck::SdfWriteError::SubstanceGroup(v)=>Some(v),_=>None}}
});

type TopologySource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    crate::error_framework::ErrorObject<TemplateAttachmentOrderError>,
>;
pub(crate) type MolWriteSource = Source<
    crate::error_framework::ErrorObject<ValenceError>,
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
                                crate::error_framework::ErrorObject<MappingValidationError>,
                                Source<
                                    crate::error_framework::ErrorObject<UIntLexicalReadError>,
                                    SdfWriteSource,
                                >,
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
#[doc = "Structured QueryGraphError payload of MolWriteError.QueryGraph; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_graph(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<QueryGraphError>>>{match self.inner.as_ref(){ck::MolWriteError::QueryGraph(value)=>Ok(Some(QueryGraphError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured QueryAtomConversionError payload of MolWriteError.QueryAtom; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_atom(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<QueryAtomConversionError>>>{match self.inner.as_ref(){ck::MolWriteError::QueryAtom(value)=>Ok(Some(QueryAtomConversionError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured QueryStateError payload of MolWriteError.QueryState; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_state(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<QueryStateError>>>{match self.inner.as_ref(){ck::MolWriteError::QueryState(value)=>Ok(Some(QueryStateError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured SmartsWriteError payload of MolWriteError.QuerySmarts; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_smarts(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<SmartsWriteError>>>{match self.inner.as_ref(){ck::MolWriteError::QuerySmarts(value)=>Ok(Some(SmartsWriteError::from_ref(self.inner.node().project(|e|match e.downcast_ref::<ck::MolWriteError>().unwrap(){ck::MolWriteError::QuerySmarts(value)=>value,_=>unreachable!()})).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured ValenceError payload of MolWriteError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>>{match self.inner.as_ref(){ck::MolWriteError::Valence(value)=>Ok(Some(ValenceError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured KekulizeError payload of MolWriteError.Kekulize; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn kekulize(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<KekulizeError>>>{match self.inner.as_ref(){ck::MolWriteError::Kekulize(value)=>Ok(Some(KekulizeError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured AtropisomerError payload of MolWriteError.Atropisomer; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn atropisomer(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<AtropisomerError>>>{match self.inner.as_ref(){ck::MolWriteError::Atropisomer(value)=>Ok(Some(AtropisomerError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured WedgeError payload of MolWriteError.Wedge; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn wedge(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<WedgeError>>>{match self.inner.as_ref(){ck::MolWriteError::Wedge(value)=>Ok(Some(WedgeError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured Coordinate2DError payload of MolWriteError.Depict; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn depict(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<Coordinate2DError>>>{match self.inner.as_ref(){ck::MolWriteError::Depict(value)=>Ok(Some(Coordinate2DError::from_ref(self.inner.node().project(|e|match e.downcast_ref::<ck::MolWriteError>().unwrap(){ck::MolWriteError::Depict(value)=>value,_=>unreachable!()})).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured TopologyValidationError payload of MolWriteError.Topology; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn topology(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>>{match self.inner.as_ref(){ck::MolWriteError::Topology(value)=>Ok(Some(TopologyValidationError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured CoordinateValidationError payload of MolWriteError.Coordinates; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn coordinates(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<CoordinateValidationError>>>{match self.inner.as_ref(){ck::MolWriteError::Coordinates(value)=>Ok(Some(CoordinateValidationError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured PropertyValueError payload of MolWriteError.Property; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn property_error(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PropertyValueError>>>{match self.inner.as_ref(){ck::MolWriteError::Property(value)=>Ok(Some(PropertyValueError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured PropertyUIntReadError payload of MolWriteError.UnsignedProperty; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn unsigned_property(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PropertyUIntReadError>>>{match self.inner.as_ref(){ck::MolWriteError::UnsignedProperty(value)=>Ok(Some(PropertyUIntReadError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
#[doc = "Structured SdfWriteError payload of MolWriteError.Detached; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn detached(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<SdfWriteError>>>{match self.inner.as_ref(){ck::MolWriteError::Detached(value)=>Ok(Some(SdfWriteError::from_rust(value.clone()).value(py).extract()?)),_=>Ok(None)}}
#[getter]#[doc = "Capability required by the failing native operation. Applies to MolWriteError.MissingCapability. None when this payload is absent."] fn capability(&self)->Option<&'static str>{match self.inner.as_ref(){ck::MolWriteError::MissingCapability(v)=>Some(*v),_=>None}}
#[getter]#[doc="Unsupported molecular output subset; None for other writer failures."]fn unsupported_subset(&self)->Option<&'static str>{match self.inner.as_ref(){ck::MolWriteError::UnsupportedSubset(v)=>Some(*v),_=>None}}
#[getter]#[doc = "Value rejected or reported by the native operation. Applies to MolWriteError.Value. None when this payload is absent."] fn value(&self)->Option<&str>{match self.inner.as_ref(){ck::MolWriteError::Value(v)=>Some(v),_=>None}}
});

type MolecularIoSource = Source<
    crate::error_framework::ErrorObject<IoError>,
    Source<
        crate::error_framework::ErrorObject<XyzReadError>,
        Source<
            crate::error_framework::ErrorObject<XyzWriteError>,
            Source<
                crate::error_framework::ErrorObject<OperationError>,
                Source<
                    crate::error_framework::ErrorObject<Mol2ReadError>,
                    Source<
                        crate::error_framework::ErrorObject<Mol2PostError>,
                        Source<
                            crate::error_framework::ErrorObject<MolWriteError>,
                            Source<
                                crate::error_framework::ErrorObject<PdbWriteError>,
                                Source<
                                    crate::error_framework::ErrorObject<BioPdbReadError>,
                                    Source<
                                        crate::error_framework::ErrorObject<BioMmcifReadError>,
                                        Source<
                                            crate::error_framework::ErrorObject<BioMoleculeError>,
                                            Source<
                                                crate::error_framework::ErrorObject<SdfError>,
                                                crate::error_framework::ErrorObject<FromUtf8Error>,
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
#[getter]#[doc = "Filesystem path involved in the operation. Applies to MolecularIoError.Io. None when this payload is absent."] fn path(&self)->Option<std::path::PathBuf>{match self.inner.as_ref(){ck::MolecularIoError::Io{path,..}=>Some(path.clone()),_=>None}}
#[getter]#[doc = "Format involved in the failed conversion. Applies to MolecularIoError.NoRecord, MolecularIoError.OutputUtf8. None when this payload is absent."] fn format(&self)->Option<&'static str>{match self.inner.as_ref(){ck::MolecularIoError::NoRecord{format}|ck::MolecularIoError::OutputUtf8{format,..}=>Some(*format),_=>None}}
});

error_projection!(retained AvalonFingerprintError,AvalonFingerprintErrorKind,ck::AvalonFingerprintError,"fingerprint",{
Input(_),Engine(_),
},cause=|inner:&RetainedError<ck::AvalonFingerprintError>|{
let node=inner.node().source()?;
Some(match inner.as_ref(){ck::AvalonFingerprintError::Input(_)=>MolecularIoError::from_ref(node),ck::AvalonFingerprintError::Engine(v)=>AvalonEngineError::from_rust(v.clone())})
},source=Option<Source<crate::error_framework::ErrorObject<MolecularIoError>,crate::error_framework::ErrorObject<AvalonEngineError>>>);
