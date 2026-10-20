//! Exact readonly mmCIF reader causes, retaining their original native owner.
#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::types::PyModule;
    use std::sync::Arc;

    #[test]
    fn mmcif_factory_and_complete_causes_keep_native_fields_bits_sources_and_ownership() {
        Python::initialize();
        Python::attach(|py|->PyResult<()> {
            let text=concat!("data_example\nloop_\n",
                "_atom_site.id\n_atom_site.group_PDB\n_atom_site.type_symbol\n",
                "_atom_site.label_atom_id\n_atom_site.label_alt_id\n_atom_site.label_comp_id\n",
                "_atom_site.label_asym_id\n_atom_site.label_entity_id\n_atom_site.label_seq_id\n",
                "_atom_site.pdbx_PDB_ins_code\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n",
                "_atom_site.occupancy\n_atom_site.B_iso_or_equiv\n_atom_site.pdbx_formal_charge\n",
                "_atom_site.auth_seq_id\n_atom_site.auth_comp_id\n_atom_site.auth_asym_id\n_atom_site.auth_atom_id\n",
                "1 ATOM C CA . GLY L1 E1 100 . -0.0 1.25 2.5 1 20 0 10 GLY A CA\n");
            let expected=ck::BioStructure::from_mmcif(text).unwrap();
            let api=PyModule::new(py,"mmcif_reader_test")?;
            api.add_function(wrap_pyfunction!(bio_from_mmcif,&api)?)?;
            let factory=py.get_type::<BioStructure>();
            for actual in [factory.call_method1("from_mmcif",(text,))?,api.call_method1("bio_from_mmcif",(text,))?] {
                let actual=actual.extract::<Py<BioStructure>>()?;
                assert_eq!(actual.borrow(py).inner,expected);
                let bits=|value:&ck::BioStructure|value.coordinates().positions().iter().map(|row|row.map(f64::to_bits)).collect::<Vec<_>>();
                assert_eq!(bits(&actual.borrow(py).inner),bits(&expected));
            }
            let invalid="data_x\n_incomplete";
            let native=ck::BioStructure::from_mmcif(invalid).unwrap_err();
            let display=native.to_string();
            let cif=native.source().unwrap().downcast_ref::<ck::CifReadError>().unwrap().clone();
            let error=factory.call_method1("from_mmcif",(invalid,)).unwrap_err();
            assert!(error.is_instance_of::<BioMmcifReadError>(py));
            assert_eq!(error.value(py).str()?.to_str()?,display);
            assert!(error.value(py).call_method0("stage")?.extract::<BioMmcifReadStage>()?==native.stage().into());
            let source=error.value(py).getattr("source")?.unbind();
            assert!(source.bind(py).is(error.value(py).getattr("__cause__")?.as_any()));
            assert!(source.bind(py).is_instance_of::<CifReadError>());
            assert_eq!(source.bind(py).getattr("source")?.extract::<String>()?,cif.source());
            drop(error);drop(native);
            assert_eq!(source.bind(py).getattr("line")?.extract::<usize>()?,cif.line());
            assert!(source.bind(py).getattr("__cause__")?.is_none());
            let missing=text.replace("_atom_site.label_atom_id\n","")
                .replace("_atom_site.auth_atom_id\n","")
                .replace("1 ATOM C CA . GLY","1 ATOM C . GLY")
                .replace("10 GLY A CA\n", "10 GLY A\n");
            let error=factory.call_method1("from_mmcif",(&missing,)).unwrap_err();
            let grouping=error.value(py).getattr("source")?;
            assert!(grouping.is_instance_of::<MmcifAtomSiteGroupingError>());
            let field=grouping.getattr("source")?;
            assert!(field.is_instance_of::<MissingAtomSiteIdentityField>());
            assert!(field.getattr("field")?.extract::<AtomSiteIdentityField>()?==AtomSiteIdentityField::AtomId);
            assert!(field.setattr("field",py.None()).is_err());

            let nested=ck::MmcifAssemblyError::Transform(ck::MmcifCrystalInfoError::Cif(cif.clone()));
            let root=RetainedError::new(nested);
            let native_child=root.as_ref().source().unwrap().downcast_ref::<ck::MmcifCrystalInfoError>().unwrap();
            let error=MmcifAssemblyError::from_ref(root.node().clone());
            let child=error.value(py).getattr("source")?.extract::<Py<MmcifCrystalInfoError>>()?;
            assert!(std::ptr::eq(child.borrow(py).inner.as_ref(),native_child));
            drop(error);drop(root);
            assert!(child.bind(py).getattr("source")?.is_instance_of::<CifReadError>());

            // Every reader-specific error projection remains exact, including
            // string-valued `source`, mixed native value widths and source-free leaves.
            let nodes:Vec<Arc<dyn Error+Send+Sync>>=vec![
                Arc::new(ck::CoordinateBlockSelectionError::LaterCoordinateBlock{block_number:usize::MAX,source:"input\0name".into()}),
                Arc::new(ck::MmcifCrystalInfoError::StagingInvariant{stage:"crystal"}),
                Arc::new(ck::MmcifTlsInfoError::NumericIdOutsideSourceDefinedRange{id:"huge".into()}),
                Arc::new(ck::MmcifSiftsUnpError::NumberOutsideU16{value:"huge".into()}),
                Arc::new(ck::MmcifEntityInfoError::Sifts(ck::MmcifSiftsUnpError::Cif(cif.clone()))),
                Arc::new(ck::MmcifSequenceIdError::InconsistentInsertionCode{sequence_id:"1A".into()}),
                Arc::new(ck::MmcifLabelAddressError::ChainNameNotRepresentable{value:"chain".into()}),
                Arc::new(ck::MmcifConnectionError::SymmetryOperatorNotRepresentable{value:i32::MIN}),
                Arc::new(ck::MmcifCisPepError::SequenceId{partner:"begin",error:ck::MmcifSequenceIdError::InvalidSequenceNumber(cif.clone())}),
                Arc::new(ck::MmcifModResError::ResidueAddressNotRepresentable),
                Arc::new(ck::MmcifOperationExprError::SourceUndefinedRangeIncrement),
                Arc::new(ck::MmcifAssemblyError::OperationExpr(ck::MmcifOperationExprError::SourceUndefinedIntegerOverflow)),
                Arc::new(ck::MmcifAtomSiteScalarError::SerialOutsideSourceDefinedRange{value:"huge".into()}),
                Arc::new(ck::MmcifAtomSiteGroupingError::UndefinedGroupPdbIndex{value:"".into(),index:usize::MAX,length:0}),
                Arc::new(ck::MmcifHelixError::SequenceId{endpoint:"start",error:ck::MmcifSequenceIdError::InvalidInsertionCode(cif.clone())}),
                Arc::new(ck::MmcifSheetError::ChainNameNotRepresentable{address:"begin",value:"chain".into()}),
                Arc::new(ck::MmcifCcdRestoreError::SequenceUtf8(String::from_utf8(vec![0xff]).unwrap_err())),
                Arc::new(ck::PdbResidueKeyError::ResidueNameNotAscii{byte:0xff}),
            ];
            for node in nodes {
                let display=node.to_string();
                let has_source=node.source().is_some();
                let error=project(ErrorRef::new(node));
                assert_eq!(error.value(py).str()?.to_str()?,display);
                assert_eq!(error.cause(py).is_some(),has_source);
                assert!(error.value(py).setattr("variant",py.None()).is_err());
            }
            let error=project(ErrorRef::new(Arc::new(ck::CoordinateBlockSelectionError::LaterCoordinateBlock{block_number:usize::MAX,source:"input\0name".into()})));
            assert_eq!(error.value(py).getattr("source")?.extract::<String>()?,"input\0name");
            assert_eq!(error.value(py).getattr("block_number")?.extract::<usize>()?,usize::MAX);
            assert!(error.cause(py).is_none());
            let error=project(ErrorRef::new(Arc::new(ck::MmcifConnectionError::SymmetryOperatorNotRepresentable{value:i32::MIN})));
            assert_eq!(error.value(py).getattr("value")?.extract::<i32>()?,i32::MIN);
            let error=project(ErrorRef::new(Arc::new(ck::MmcifHelixError::SequenceId{endpoint:"start",error:ck::MmcifSequenceIdError::InvalidInsertionCode(cif)})));
            assert!(error.value(py).getattr("error")?.is(&error.value(py).getattr("__cause__")?));
            Ok(())
        }).unwrap();
    }
}
use crate::bio_errors::BioStructureError;
use crate::bio_mmcif_writer::CifReadError;
use crate::bio_structure::BioStructure;
use crate::chemical_values::enum_value;
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::ErrorRef;
use crate::public_selectors::BioMmcifReadStage;
use crate::retained_error::RetainedError;
use crate::utf8_errors::FromUtf8Error;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

fn project(node: ErrorRef) -> PyErr {
    macro_rules! retained {($($name:ident),+)=>{$(if node.error().is::<ck::$name>() {return $name::from_ref(node);})+};}
    retained!(
        CoordinateBlockSelectionError,
        MmcifCrystalInfoError,
        MmcifTlsInfoError,
        MmcifHelixError,
        MmcifSheetError,
        MmcifSiftsUnpError,
        MmcifEntityInfoError,
        MmcifSequenceIdError,
        MmcifAtomSiteGroupingError,
        MissingAtomSiteIdentityField,
        MmcifLabelAddressError,
        MmcifConnectionError,
        MmcifCisPepError,
        MmcifModResError,
        MmcifOperationExprError,
        MmcifAssemblyError,
        MmcifAtomSiteScalarError,
        MmcifAtomSiteRowOutOfBounds,
        MmcifCcdRestoreError,
        PdbResidueKeyError
    );
    if let Some(value) = node.error().downcast_ref::<ck::CifReadError>() {
        return CifReadError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::BioStructureError>() {
        return BioStructureError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<std::string::FromUtf8Error>() {
        return FromUtf8Error::from_rust(value.clone());
    }
    unreachable!("unmodeled native mmCIF reader cause")
}

error_projection!(retained PdbResidueKeyError,PdbResidueKeyErrorKind,ck::PdbResidueKeyError,"bio",{
    SequenceIdNotRepresentable,ResidueNameTooWide{width},ResidueNameNotAscii{byte},
    SegmentTooWide{width},SegmentNotAscii{byte},SegmentNotRepresentable{width},
},cause=|_:&RetainedError<ck::PdbResidueKeyError>|None::<PyErr>,methods={
    #[getter] fn width(&self)->Option<usize>{match self.inner.as_ref(){ck::PdbResidueKeyError::ResidueNameTooWide{width}|ck::PdbResidueKeyError::SegmentTooWide{width}|ck::PdbResidueKeyError::SegmentNotRepresentable{width}=>Some(*width),_=>None}}
    #[getter] fn byte(&self)->Option<u8>{match self.inner.as_ref(){ck::PdbResidueKeyError::ResidueNameNotAscii{byte}|ck::PdbResidueKeyError::SegmentNotAscii{byte}=>Some(*byte),_=>None}}
});
// `source` is the actual input-name string here, not an Error::source edge.
error_projection!(retained CoordinateBlockSelectionError,CoordinateBlockSelectionErrorKind,ck::CoordinateBlockSelectionError,"bio",{
    EmptyDocument{source},LaterCoordinateBlock{block_number:usize,source},
},cause=|_:&RetainedError<ck::CoordinateBlockSelectionError>|None::<PyErr>,methods={
    #[getter] fn source(&self)->&str{match self.inner.as_ref(){ck::CoordinateBlockSelectionError::EmptyDocument{source}|ck::CoordinateBlockSelectionError::LaterCoordinateBlock{source,..}=>source}}
});
error_projection!(retained MmcifCrystalInfoError,MmcifCrystalInfoErrorKind,ck::MmcifCrystalInfoError,"bio",{
    Cif(_),Bio(_),StagingInvariant{stage:&'static str},
},cause=|inner:&RetainedError<ck::MmcifCrystalInfoError>|inner.node().source().map(project),source=Option<Source<Py<CifReadError>,Py<BioStructureError>>>);
error_projection!(retained MmcifSequenceIdError,MmcifSequenceIdErrorKind,ck::MmcifSequenceIdError,"bio",{
    InconsistentInsertionCode{sequence_id:String},InvalidInsertionCode(_),InvalidSequenceNumber(_),
},cause=|inner:&RetainedError<ck::MmcifSequenceIdError>|inner.node().source().map(project),source=Option<Py<CifReadError>>);
error_projection!(retained MmcifTlsInfoError,MmcifTlsInfoErrorKind,ck::MmcifTlsInfoError,"bio",{
    Cif(_),SequenceId(_),NumericIdOutsideSourceDefinedRange{id:String},SelectionChainNotRepresentable{chain:String},
},cause=|inner:&RetainedError<ck::MmcifTlsInfoError>|inner.node().source().map(project),source=Option<Source<Py<CifReadError>,Py<MmcifSequenceIdError>>>);
error_projection!(retained MmcifSiftsUnpError,MmcifSiftsUnpErrorKind,ck::MmcifSiftsUnpError,"bio",{
    Cif(_),EntityNotFound{entity_id:String},SubchainNotFound{asym_id:String},SequenceNotFound{seq_id:String},NumberOutsideU16{value:String},StagingInvariant{stage:&'static str},
},cause=|inner:&RetainedError<ck::MmcifSiftsUnpError>|inner.node().source().map(project),source=Option<Py<CifReadError>>);
error_projection!(retained MmcifEntityInfoError,MmcifEntityInfoErrorKind,ck::MmcifEntityInfoError,"bio",{
    Cif(_),SequenceId(_),Sifts(_),
},cause=|inner:&RetainedError<ck::MmcifEntityInfoError>|inner.node().source().map(project),source=Option<Source<Py<CifReadError>,Source<Py<MmcifSequenceIdError>,Py<MmcifSiftsUnpError>>>>);
error_projection!(retained MmcifLabelAddressError,MmcifLabelAddressErrorKind,ck::MmcifLabelAddressError,"bio",{
    LabelSequenceId(_),ChainNameNotRepresentable{value:String},
},cause=|inner:&RetainedError<ck::MmcifLabelAddressError>|inner.node().source().map(project),source=Option<Py<CifReadError>>);
error_projection!(retained MmcifModResError,MmcifModResErrorKind,ck::MmcifModResError,"bio",{
    Cif(_),ChainNameNotRepresentable{value},ResidueNameNotRepresentable{value},ResidueAddressNotRepresentable,SequenceId(_),
},cause=|inner:&RetainedError<ck::MmcifModResError>|inner.node().source().map(project),source=Option<Source<Py<CifReadError>,Py<MmcifSequenceIdError>>>,methods={
    #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){ck::MmcifModResError::ChainNameNotRepresentable{value}|ck::MmcifModResError::ResidueNameNotRepresentable{value}=>Some(value),_=>None}}
});
error_projection!(retained MmcifOperationExprError,MmcifOperationExprErrorKind,ck::MmcifOperationExprError,"bio",{
    SourceUndefinedIntegerOverflow,SourceUndefinedRangeIncrement,
},cause=|_:&RetainedError<ck::MmcifOperationExprError>|None::<PyErr>);
error_projection!(retained MmcifAssemblyError,MmcifAssemblyErrorKind,ck::MmcifAssemblyError,"bio",{
    Cif(_),Transform(_),OperationExpr(_),
},cause=|inner:&RetainedError<ck::MmcifAssemblyError>|inner.node().source().map(project),source=Option<Source<Py<CifReadError>,Source<Py<MmcifCrystalInfoError>,Py<MmcifOperationExprError>>>>);
error_projection!(retained MmcifAtomSiteScalarError,MmcifAtomSiteScalarErrorKind,ck::MmcifAtomSiteScalarError,"bio",{
    SerialOutsideSourceDefinedRange{value},TlsOutsideSourceDefinedRange{value},CifCharacter(_),CifInteger(_),CifNumber(_),
},cause=|inner:&RetainedError<ck::MmcifAtomSiteScalarError>|inner.node().source().map(project),source=Option<Py<CifReadError>>,methods={
    #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){ck::MmcifAtomSiteScalarError::SerialOutsideSourceDefinedRange{value}|ck::MmcifAtomSiteScalarError::TlsOutsideSourceDefinedRange{value}=>Some(value),_=>None}}
});
error_projection!(retained MmcifCcdRestoreError,MmcifCcdRestoreErrorKind,ck::MmcifCcdRestoreError,"bio",{
    Residue(_),SequenceUtf8(_),
},cause=|inner:&RetainedError<ck::MmcifCcdRestoreError>|inner.node().source().map(project),source=Option<Source<Py<PdbResidueKeyError>,Py<FromUtf8Error>>>);

error_projection!(retained MmcifHelixError,MmcifHelixErrorKind,ck::MmcifHelixError,"bio",{
 Cif(_),ChainNameNotRepresentable{endpoint,value},ResidueNameNotRepresentable{endpoint,value},
 ResidueAddressNotRepresentable{endpoint},SequenceId{endpoint,error},HelixClass(_),HelixLength(_),
},cause=|inner:&RetainedError<ck::MmcifHelixError>|inner.node().source().map(project),
source=Option<Source<Py<CifReadError>,Py<MmcifSequenceIdError>>>,methods={
 #[getter] fn endpoint(&self)->Option<&'static str>{match self.inner.as_ref(){
 ck::MmcifHelixError::ChainNameNotRepresentable{endpoint,..}|ck::MmcifHelixError::ResidueNameNotRepresentable{endpoint,..}
 |ck::MmcifHelixError::ResidueAddressNotRepresentable{endpoint}|ck::MmcifHelixError::SequenceId{endpoint,..}=>Some(*endpoint),_=>None}}
 #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){
 ck::MmcifHelixError::ChainNameNotRepresentable{value,..}|ck::MmcifHelixError::ResidueNameNotRepresentable{value,..}=>Some(value),_=>None}}
 #[getter] #[doc=r#"{"role":"error_payload"}"#]
 fn error(slf:Bound<'_,Self>)->PyResult<Option<Py<MmcifSequenceIdError>>>{
 if matches!(slf.borrow().inner.as_ref(),ck::MmcifHelixError::SequenceId{..}) {Ok(Some(slf.getattr("__cause__")?.extract()?))}else{Ok(None)}
 }
});

error_projection!(retained MmcifSheetError,MmcifSheetErrorKind,ck::MmcifSheetError,"bio",{
 Cif(_),ChainNameNotRepresentable{address,value},ResidueNameNotRepresentable{address,value},
 ResidueAddressNotRepresentable{address},SequenceId{address,error},
},cause=|inner:&RetainedError<ck::MmcifSheetError>|inner.node().source().map(project),
source=Option<Source<Py<CifReadError>,Py<MmcifSequenceIdError>>>,methods={
 #[getter] fn address(&self)->Option<&'static str>{match self.inner.as_ref(){
 ck::MmcifSheetError::ChainNameNotRepresentable{address,..}|ck::MmcifSheetError::ResidueNameNotRepresentable{address,..}
 |ck::MmcifSheetError::ResidueAddressNotRepresentable{address}|ck::MmcifSheetError::SequenceId{address,..}=>Some(*address),_=>None}}
 #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){
 ck::MmcifSheetError::ChainNameNotRepresentable{value,..}|ck::MmcifSheetError::ResidueNameNotRepresentable{value,..}=>Some(value),_=>None}}
 #[getter] #[doc=r#"{"role":"error_payload"}"#]
 fn error(slf:Bound<'_,Self>)->PyResult<Option<Py<MmcifSequenceIdError>>>{
 if matches!(slf.borrow().inner.as_ref(),ck::MmcifSheetError::SequenceId{..}) {Ok(Some(slf.getattr("__cause__")?.extract()?))}else{Ok(None)}
 }
});

error_projection!(retained MmcifCisPepError,MmcifCisPepErrorKind,ck::MmcifCisPepError,"bio",{
 Cif(_),ChainNameNotRepresentable{partner,value},ResidueNameNotRepresentable{partner,value},
 ResidueAddressNotRepresentable{partner},SequenceId{partner,error},
},cause=|inner:&RetainedError<ck::MmcifCisPepError>|inner.node().source().map(project),
source=Option<Source<Py<CifReadError>,Py<MmcifSequenceIdError>>>,methods={
 #[getter] fn partner(&self)->Option<&'static str>{match self.inner.as_ref(){
 ck::MmcifCisPepError::ChainNameNotRepresentable{partner,..}|ck::MmcifCisPepError::ResidueNameNotRepresentable{partner,..}
 |ck::MmcifCisPepError::ResidueAddressNotRepresentable{partner}|ck::MmcifCisPepError::SequenceId{partner,..}=>Some(*partner),_=>None}}
 #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){
 ck::MmcifCisPepError::ChainNameNotRepresentable{value,..}|ck::MmcifCisPepError::ResidueNameNotRepresentable{value,..}=>Some(value),_=>None}}
 #[getter] #[doc=r#"{"role":"error_payload"}"#]
 fn error(slf:Bound<'_,Self>)->PyResult<Option<Py<MmcifSequenceIdError>>>{
 if matches!(slf.borrow().inner.as_ref(),ck::MmcifCisPepError::SequenceId{..}) {Ok(Some(slf.getattr("__cause__")?.extract()?))}else{Ok(None)}
 }
});

enum_value!(AtomSiteIdentityField,{AsymId,CompId,AtomId,SeqId});

error_projection!(retained MmcifConnectionError,MmcifConnectionErrorKind,ck::MmcifConnectionError,"bio",{
 Cif(_),SequenceId(_),LabelAddress(_),MissingAddressIdentifiers{partner},NoStructuralModels,
 ChainNameNotRepresentable{partner,value},ResidueNameNotRepresentable{partner,value},ResidueAddressNotRepresentable{partner,value},
 SymmetryOperatorOutsideSourceDefinedRange{value},SymmetryOperatorNotRepresentable{value},
},cause=|inner:&RetainedError<ck::MmcifConnectionError>|inner.node().source().map(project),source=Option<Source<Py<CifReadError>,Source<Py<MmcifSequenceIdError>,Py<MmcifLabelAddressError>>>>,methods={
 #[getter] fn partner(&self)->Option<usize>{match self.inner.as_ref(){
 ck::MmcifConnectionError::MissingAddressIdentifiers{partner}|ck::MmcifConnectionError::ChainNameNotRepresentable{partner,..}
 |ck::MmcifConnectionError::ResidueNameNotRepresentable{partner,..}|ck::MmcifConnectionError::ResidueAddressNotRepresentable{partner,..}=>Some(*partner),_=>None}}
 #[getter] fn value(&self)->Option<Source<String,i32>>{match self.inner.as_ref(){
 ck::MmcifConnectionError::ChainNameNotRepresentable{value,..}|ck::MmcifConnectionError::ResidueNameNotRepresentable{value,..}
 |ck::MmcifConnectionError::ResidueAddressNotRepresentable{value,..}|ck::MmcifConnectionError::SymmetryOperatorOutsideSourceDefinedRange{value}=>Some(Source::First(value.clone())),
 ck::MmcifConnectionError::SymmetryOperatorNotRepresentable{value}=>Some(Source::Second(*value)),_=>None}}
});
type GroupingSource = Source<
    Py<MissingAtomSiteIdentityField>,
    Source<Py<MmcifSequenceIdError>, Source<Py<MmcifAtomSiteScalarError>, Py<CifReadError>>>,
>;
error_projection!(retained MmcifAtomSiteGroupingError,MmcifAtomSiteGroupingErrorKind,ck::MmcifAtomSiteGroupingError,"bio",{
 IdentityField(_),SequenceId(_),AtomNameNotRepresentable{value},AtomScalar(_),ModelNumber(_),LabelSequenceNumber(_),
 AuthorChainIdNotRepresentable{value},ResidueNameNotRepresentable{value},ResidueAddressNotRepresentable{value},ResidueSourceIdsNotRepresentable,
 UndefinedGroupPdbIndex{value,index:usize,length:usize},MissingActiveModel{row_index:usize},
},cause=|inner:&RetainedError<ck::MmcifAtomSiteGroupingError>|inner.node().source().map(project),source=Option<GroupingSource>,methods={
 #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){
 ck::MmcifAtomSiteGroupingError::AtomNameNotRepresentable{value}|ck::MmcifAtomSiteGroupingError::AuthorChainIdNotRepresentable{value}
 |ck::MmcifAtomSiteGroupingError::ResidueNameNotRepresentable{value}|ck::MmcifAtomSiteGroupingError::ResidueAddressNotRepresentable{value}
 |ck::MmcifAtomSiteGroupingError::UndefinedGroupPdbIndex{value,..}=>Some(value),_=>None}}
});

type MmcifReadSource = Source<
    Py<CifReadError>,
    Source<
        Py<BioStructureError>,
        Source<
            Py<CoordinateBlockSelectionError>,
            Source<
                Py<MmcifCrystalInfoError>,
                Source<
                    Py<MmcifTlsInfoError>,
                    Source<
                        Py<MmcifAtomSiteGroupingError>,
                        Source<
                            Py<MmcifEntityInfoError>,
                            Source<
                                Py<MmcifHelixError>,
                                Source<
                                    Py<MmcifSheetError>,
                                    Source<
                                        Py<MmcifConnectionError>,
                                        Source<
                                            Py<MmcifCisPepError>,
                                            Source<
                                                Py<MmcifModResError>,
                                                Source<
                                                    Py<MmcifAssemblyError>,
                                                    Source<
                                                        Py<MmcifSiftsUnpError>,
                                                        Source<
                                                            Py<MmcifCcdRestoreError>,
                                                            Py<MmcifAtomSiteRowOutOfBounds>,
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
            >,
        >,
    >,
>;
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct BioMmcifReadError {
    inner: RetainedError<ck::BioMmcifReadError>,
}
impl BioMmcifReadError {
    pub(crate) fn from_rust(inner: ck::BioMmcifReadError) -> PyErr {
        Self::from_stored(RetainedError::new(inner))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Self::from_stored(RetainedError::from_ref(node))
    }
    fn from_stored(inner: RetainedError<ck::BioMmcifReadError>) -> PyErr {
        Python::attach(|py| {
            let cause = inner.node().source().map(project);
            match Py::new(py, Self { inner }) {
                Ok(value) => {
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl BioMmcifReadError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn stage(&self) -> BioMmcifReadStage {
        self.inner.stage().into()
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<MmcifReadSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[pyfunction]
pub fn bio_from_mmcif(py: Python<'_>, text: &str) -> PyResult<BioStructure> {
    py.detach(|| ck::BioStructure::from_mmcif(text))
        .map(|inner| BioStructure { inner })
        .map_err(BioMmcifReadError::from_rust)
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct MissingAtomSiteIdentityField {
    inner: RetainedError<ck::MissingAtomSiteIdentityField>,
}
impl MissingAtomSiteIdentityField {
    fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl MissingAtomSiteIdentityField {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn field(&self) -> AtomSiteIdentityField {
        self.inner.field().into()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct MmcifAtomSiteRowOutOfBounds {
    inner: RetainedError<ck::MmcifAtomSiteRowOutOfBounds>,
}
impl MmcifAtomSiteRowOutOfBounds {
    fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl MmcifAtomSiteRowOutOfBounds {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn row_index(&self) -> usize {
        self.inner.row_index()
    }
    #[getter]
    fn parsed_count(&self) -> usize {
        self.inner.parsed_count()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
