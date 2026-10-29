//! Native BIO/MMCIF errors retain their original owner and payloads.
use crate::errors::{enum_error, error_projection, exception};
use crate::external_errors::FromUtf8Error;
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
enum_error!(PdbResidueKeyError,PdbResidueKeyErrorKind,ck::PdbResidueKeyError,{
SequenceIdNotRepresentable=>ck::PdbResidueKeyError::SequenceIdNotRepresentable,
ResidueNameTooWide=>ck::PdbResidueKeyError::ResidueNameTooWide{..},
ResidueNameNotAscii=>ck::PdbResidueKeyError::ResidueNameNotAscii{..},
SegmentTooWide=>ck::PdbResidueKeyError::SegmentTooWide{..},
SegmentNotAscii=>ck::PdbResidueKeyError::SegmentNotAscii{..},
SegmentNotRepresentable=>ck::PdbResidueKeyError::SegmentNotRepresentable{..}
},{
#[wasm_bindgen(getter,js_name=width)] pub fn width(&self)->Option<usize>{match self.inner.as_ref(){
ck::PdbResidueKeyError::ResidueNameTooWide{width,..}=>Some(*width),
ck::PdbResidueKeyError::SegmentTooWide{width,..}=>Some(*width),
ck::PdbResidueKeyError::SegmentNotRepresentable{width,..}=>Some(*width),
_=>None,}}
#[wasm_bindgen(getter,js_name=byte)] pub fn byte(&self)->Option<u8>{match self.inner.as_ref(){
ck::PdbResidueKeyError::ResidueNameNotAscii{byte,..}=>Some(*byte),
ck::PdbResidueKeyError::SegmentNotAscii{byte,..}=>Some(*byte),
_=>None,}}
});
enum_error!(CoordinateBlockSelectionError,CoordinateBlockSelectionErrorKind,ck::CoordinateBlockSelectionError,{
EmptyDocument=>ck::CoordinateBlockSelectionError::EmptyDocument{..},
LaterCoordinateBlock=>ck::CoordinateBlockSelectionError::LaterCoordinateBlock{..}
},{
#[wasm_bindgen(getter,js_name=source)] pub fn source(&self)->Option<String>{match self.inner.as_ref(){
ck::CoordinateBlockSelectionError::EmptyDocument{source,..}=>Some(source.to_string()),
ck::CoordinateBlockSelectionError::LaterCoordinateBlock{source,..}=>Some(source.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=blockNumber)] pub fn block_number(&self)->Option<usize>{match self.inner.as_ref(){
ck::CoordinateBlockSelectionError::LaterCoordinateBlock{block_number,..}=>Some(*block_number),
_=>None,}}
});
enum_error!(BioStructureError,BioStructureErrorKind,ck::BioStructureError,{
RowIndexTooLarge=>ck::BioStructureError::RowIndexTooLarge{..},
RowSpanOverflow=>ck::BioStructureError::RowSpanOverflow{..},
RowSpanOutOfBounds=>ck::BioStructureError::RowSpanOutOfBounds{..},
TableTooLarge=>ck::BioStructureError::TableTooLarge{..},
NonContiguousSpan=>ck::BioStructureError::NonContiguousSpan{..},
IncompleteCoverage=>ck::BioStructureError::IncompleteCoverage{..},
ParentMismatch=>ck::BioStructureError::ParentMismatch{..},
RowReferenceOutOfBounds=>ck::BioStructureError::RowReferenceOutOfBounds{..},
CoordinateCountMismatch=>ck::BioStructureError::CoordinateCountMismatch{..},
EntitySubchainMismatch=>ck::BioStructureError::EntitySubchainMismatch{..},
EmptyResidueSpan=>ck::BioStructureError::EmptyResidueSpan{..},
ImpossibleCrystalAngle=>ck::BioStructureError::ImpossibleCrystalAngle,
AtomNotFound=>ck::BioStructureError::AtomNotFound
},{
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::RowIndexTooLarge{value,..}=>Some(*value),
_=>None,}}
#[wasm_bindgen(getter,js_name=start)] pub fn start(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::RowSpanOverflow{start,..}=>Some(*start as usize),
ck::BioStructureError::RowSpanOutOfBounds{start,..}=>Some(*start as usize),
_=>None,}}
#[wasm_bindgen(getter,js_name=len)] pub fn len(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::RowSpanOverflow{len,..}=>Some(*len as usize),
ck::BioStructureError::RowSpanOutOfBounds{len,..}=>Some(*len as usize),
ck::BioStructureError::TableTooLarge{len,..}=>Some(*len),
_=>None,}}
#[wasm_bindgen(getter,js_name=tableLen)] pub fn table_len(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::RowSpanOutOfBounds{table_len,..}=>Some(*table_len),
ck::BioStructureError::IncompleteCoverage{table_len,..}=>Some(*table_len),
ck::BioStructureError::RowReferenceOutOfBounds{table_len,..}=>Some(*table_len),
_=>None,}}
#[wasm_bindgen(getter,js_name=table)] pub fn table(&self)->Option<String>{match self.inner.as_ref(){
ck::BioStructureError::TableTooLarge{table,..}=>Some(table.to_string()),
ck::BioStructureError::NonContiguousSpan{table,..}=>Some(table.to_string()),
ck::BioStructureError::IncompleteCoverage{table,..}=>Some(table.to_string()),
ck::BioStructureError::ParentMismatch{table,..}=>Some(table.to_string()),
ck::BioStructureError::RowReferenceOutOfBounds{table,..}=>Some(table.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=expectedStart)] pub fn expected_start(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::NonContiguousSpan{expected_start,..}=>Some(*expected_start as usize),
_=>None,}}
#[wasm_bindgen(getter,js_name=actualStart)] pub fn actual_start(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::NonContiguousSpan{actual_start,..}=>Some(*actual_start as usize),
_=>None,}}
#[wasm_bindgen(getter,js_name=covered)] pub fn covered(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::IncompleteCoverage{covered,..}=>Some(*covered as usize),
_=>None,}}
#[wasm_bindgen(getter,js_name=index)] pub fn index(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::ParentMismatch{index,..}=>Some(*index as usize),
ck::BioStructureError::RowReferenceOutOfBounds{index,..}=>Some(*index as usize),
_=>None,}}
#[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::CoordinateCountMismatch{atom_count,..}=>Some(*atom_count),
_=>None,}}
#[wasm_bindgen(getter,js_name=coordinateCount)] pub fn coordinate_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::CoordinateCountMismatch{coordinate_count,..}=>Some(*coordinate_count),
_=>None,}}
#[wasm_bindgen(getter,js_name=entityId)] pub fn entity_id(&self)->Option<usize>{match self.inner.as_ref(){
ck::BioStructureError::EntitySubchainMismatch{entity_id,..}=>Some(entity_id.index()),
_=>None,}}
#[wasm_bindgen(getter,js_name=subchain)] pub fn subchain(&self)->Option<String>{match self.inner.as_ref(){
ck::BioStructureError::EntitySubchainMismatch{subchain,..}=>Some(subchain.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=operation)] pub fn operation(&self)->Option<String>{match self.inner.as_ref(){
ck::BioStructureError::EmptyResidueSpan{operation,..}=>Some(operation.to_string()),
_=>None,}}
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | BioStructureError")]
    #[derive(Clone)]
    pub type MmcifCrystalInfoErrorCause;
}
enum_error!(MmcifCrystalInfoError,MmcifCrystalInfoErrorKind,ck::MmcifCrystalInfoError,{
Cif=>ck::MmcifCrystalInfoError::Cif(_),
Bio=>ck::MmcifCrystalInfoError::Bio(_),
StagingInvariant=>ck::MmcifCrystalInfoError::StagingInvariant{..}
},{
#[wasm_bindgen(getter,js_name=stage)] pub fn stage(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifCrystalInfoError::StagingInvariant{stage,..}=>Some(stage.to_string()),
_=>None,}}
},source=MmcifCrystalInfoErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifCrystalInfoError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError")]
    #[derive(Clone)]
    pub type MmcifSequenceIdErrorCause;
}
enum_error!(MmcifSequenceIdError,MmcifSequenceIdErrorKind,ck::MmcifSequenceIdError,{
InconsistentInsertionCode=>ck::MmcifSequenceIdError::InconsistentInsertionCode{..},
InvalidInsertionCode=>ck::MmcifSequenceIdError::InvalidInsertionCode(_),
InvalidSequenceNumber=>ck::MmcifSequenceIdError::InvalidSequenceNumber(_)
},{
#[wasm_bindgen(getter,js_name=sequenceId)] pub fn sequence_id(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSequenceIdError::InconsistentInsertionCode{sequence_id,..}=>Some(sequence_id.to_string()),
_=>None,}}
},source=MmcifSequenceIdErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifSequenceIdError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | MmcifSequenceIdError")]
    #[derive(Clone)]
    pub type MmcifTlsInfoErrorCause;
}
enum_error!(MmcifTlsInfoError,MmcifTlsInfoErrorKind,ck::MmcifTlsInfoError,{
Cif=>ck::MmcifTlsInfoError::Cif(_),
SequenceId=>ck::MmcifTlsInfoError::SequenceId(_),
NumericIdOutsideSourceDefinedRange=>ck::MmcifTlsInfoError::NumericIdOutsideSourceDefinedRange{..},
SelectionChainNotRepresentable=>ck::MmcifTlsInfoError::SelectionChainNotRepresentable{..}
},{
#[wasm_bindgen(getter,js_name=id)] pub fn id(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifTlsInfoError::NumericIdOutsideSourceDefinedRange{id,..}=>Some(id.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=chain)] pub fn chain(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifTlsInfoError::SelectionChainNotRepresentable{chain,..}=>Some(chain.to_string()),
_=>None,}}
},source=MmcifTlsInfoErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifTlsInfoError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError")]
    #[derive(Clone)]
    pub type MmcifSiftsUnpErrorCause;
}
enum_error!(MmcifSiftsUnpError,MmcifSiftsUnpErrorKind,ck::MmcifSiftsUnpError,{
Cif=>ck::MmcifSiftsUnpError::Cif(_),
EntityNotFound=>ck::MmcifSiftsUnpError::EntityNotFound{..},
SubchainNotFound=>ck::MmcifSiftsUnpError::SubchainNotFound{..},
SequenceNotFound=>ck::MmcifSiftsUnpError::SequenceNotFound{..},
NumberOutsideU16=>ck::MmcifSiftsUnpError::NumberOutsideU16{..},
StagingInvariant=>ck::MmcifSiftsUnpError::StagingInvariant{..}
},{
#[wasm_bindgen(getter,js_name=entityId)] pub fn entity_id(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSiftsUnpError::EntityNotFound{entity_id,..}=>Some(entity_id.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=asymId)] pub fn asym_id(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSiftsUnpError::SubchainNotFound{asym_id,..}=>Some(asym_id.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=seqId)] pub fn seq_id(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSiftsUnpError::SequenceNotFound{seq_id,..}=>Some(seq_id.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSiftsUnpError::NumberOutsideU16{value,..}=>Some(value.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=stage)] pub fn stage(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSiftsUnpError::StagingInvariant{stage,..}=>Some(stage.to_string()),
_=>None,}}
},source=MmcifSiftsUnpErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifSiftsUnpError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | MmcifSequenceIdError | MmcifSiftsUnpError")]
    #[derive(Clone)]
    pub type MmcifEntityInfoErrorCause;
}
enum_error!(MmcifEntityInfoError,MmcifEntityInfoErrorKind,ck::MmcifEntityInfoError,{
Cif=>ck::MmcifEntityInfoError::Cif(_),
SequenceId=>ck::MmcifEntityInfoError::SequenceId(_),
Sifts=>ck::MmcifEntityInfoError::Sifts(_)
},{
},source=MmcifEntityInfoErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifEntityInfoError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError")]
    #[derive(Clone)]
    pub type MmcifLabelAddressErrorCause;
}
enum_error!(MmcifLabelAddressError,MmcifLabelAddressErrorKind,ck::MmcifLabelAddressError,{
LabelSequenceId=>ck::MmcifLabelAddressError::LabelSequenceId(_),
ChainNameNotRepresentable=>ck::MmcifLabelAddressError::ChainNameNotRepresentable{..}
},{
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifLabelAddressError::ChainNameNotRepresentable{value,..}=>Some(value.to_string()),
_=>None,}}
},source=MmcifLabelAddressErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifLabelAddressError>|native_source(owner).map(JsCast::unchecked_into));
enum_error!(MmcifOperationExprError,MmcifOperationExprErrorKind,ck::MmcifOperationExprError,{
SourceUndefinedIntegerOverflow=>ck::MmcifOperationExprError::SourceUndefinedIntegerOverflow,
SourceUndefinedRangeIncrement=>ck::MmcifOperationExprError::SourceUndefinedRangeIncrement
},{
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "CifReadError | MmcifCrystalInfoError | MmcifOperationExprError"
    )]
    #[derive(Clone)]
    pub type MmcifAssemblyErrorCause;
}
enum_error!(MmcifAssemblyError,MmcifAssemblyErrorKind,ck::MmcifAssemblyError,{
Cif=>ck::MmcifAssemblyError::Cif(_),
Transform=>ck::MmcifAssemblyError::Transform(_),
OperationExpr=>ck::MmcifAssemblyError::OperationExpr(_)
},{
},source=MmcifAssemblyErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifAssemblyError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError")]
    #[derive(Clone)]
    pub type MmcifAtomSiteScalarErrorCause;
}
enum_error!(MmcifAtomSiteScalarError,MmcifAtomSiteScalarErrorKind,ck::MmcifAtomSiteScalarError,{
SerialOutsideSourceDefinedRange=>ck::MmcifAtomSiteScalarError::SerialOutsideSourceDefinedRange{..},
TlsOutsideSourceDefinedRange=>ck::MmcifAtomSiteScalarError::TlsOutsideSourceDefinedRange{..},
CifCharacter=>ck::MmcifAtomSiteScalarError::CifCharacter(_),
CifInteger=>ck::MmcifAtomSiteScalarError::CifInteger(_),
CifNumber=>ck::MmcifAtomSiteScalarError::CifNumber(_)
},{
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifAtomSiteScalarError::SerialOutsideSourceDefinedRange{value,..}=>Some(value.to_string()),
ck::MmcifAtomSiteScalarError::TlsOutsideSourceDefinedRange{value,..}=>Some(value.to_string()),
_=>None,}}
},source=MmcifAtomSiteScalarErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifAtomSiteScalarError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | MmcifSequenceIdError")]
    #[derive(Clone)]
    pub type MmcifModResErrorCause;
}
enum_error!(MmcifModResError,MmcifModResErrorKind,ck::MmcifModResError,{
Cif=>ck::MmcifModResError::Cif(_),
ChainNameNotRepresentable=>ck::MmcifModResError::ChainNameNotRepresentable{..},
ResidueNameNotRepresentable=>ck::MmcifModResError::ResidueNameNotRepresentable{..},
ResidueAddressNotRepresentable=>ck::MmcifModResError::ResidueAddressNotRepresentable,
SequenceId=>ck::MmcifModResError::SequenceId(_)
},{
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifModResError::ChainNameNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifModResError::ResidueNameNotRepresentable{value,..}=>Some(value.to_string()),
_=>None,}}
},source=MmcifModResErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifModResError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "PdbResidueKeyError | FromUtf8Error")]
    #[derive(Clone)]
    pub type MmcifCcdRestoreErrorCause;
}
enum_error!(MmcifCcdRestoreError,MmcifCcdRestoreErrorKind,ck::MmcifCcdRestoreError,{
Residue=>ck::MmcifCcdRestoreError::Residue(_),
SequenceUtf8=>ck::MmcifCcdRestoreError::SequenceUtf8(_)
},{
},source=MmcifCcdRestoreErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifCcdRestoreError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "MissingAtomSiteIdentityField | MmcifSequenceIdError | MmcifAtomSiteScalarError | CifReadError"
    )]
    #[derive(Clone)]
    pub type MmcifAtomSiteGroupingErrorCause;
}
enum_error!(MmcifAtomSiteGroupingError,MmcifAtomSiteGroupingErrorKind,ck::MmcifAtomSiteGroupingError,{
IdentityField=>ck::MmcifAtomSiteGroupingError::IdentityField(_),
SequenceId=>ck::MmcifAtomSiteGroupingError::SequenceId(_),
AtomNameNotRepresentable=>ck::MmcifAtomSiteGroupingError::AtomNameNotRepresentable{..},
AtomScalar=>ck::MmcifAtomSiteGroupingError::AtomScalar(_),
ModelNumber=>ck::MmcifAtomSiteGroupingError::ModelNumber(_),
LabelSequenceNumber=>ck::MmcifAtomSiteGroupingError::LabelSequenceNumber(_),
AuthorChainIdNotRepresentable=>ck::MmcifAtomSiteGroupingError::AuthorChainIdNotRepresentable{..},
ResidueNameNotRepresentable=>ck::MmcifAtomSiteGroupingError::ResidueNameNotRepresentable{..},
ResidueAddressNotRepresentable=>ck::MmcifAtomSiteGroupingError::ResidueAddressNotRepresentable{..},
ResidueSourceIdsNotRepresentable=>ck::MmcifAtomSiteGroupingError::ResidueSourceIdsNotRepresentable,
UndefinedGroupPdbIndex=>ck::MmcifAtomSiteGroupingError::UndefinedGroupPdbIndex{..},
MissingActiveModel=>ck::MmcifAtomSiteGroupingError::MissingActiveModel{..}
},{
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifAtomSiteGroupingError::AtomNameNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifAtomSiteGroupingError::AuthorChainIdNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifAtomSiteGroupingError::ResidueNameNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifAtomSiteGroupingError::ResidueAddressNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifAtomSiteGroupingError::UndefinedGroupPdbIndex{value,..}=>Some(value.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=index)] pub fn index(&self)->Option<usize>{match self.inner.as_ref(){
ck::MmcifAtomSiteGroupingError::UndefinedGroupPdbIndex{index,..}=>Some(*index),
_=>None,}}
#[wasm_bindgen(getter,js_name=length)] pub fn length(&self)->Option<usize>{match self.inner.as_ref(){
ck::MmcifAtomSiteGroupingError::UndefinedGroupPdbIndex{length,..}=>Some(*length),
_=>None,}}
#[wasm_bindgen(getter,js_name=rowIndex)] pub fn row_index(&self)->Option<usize>{match self.inner.as_ref(){
ck::MmcifAtomSiteGroupingError::MissingActiveModel{row_index,..}=>Some(*row_index),
_=>None,}}
},source=MmcifAtomSiteGroupingErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifAtomSiteGroupingError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "CifReadError | MmcifSequenceIdError | MmcifLabelAddressError"
    )]
    #[derive(Clone)]
    pub type MmcifConnectionErrorCause;
}
enum_error!(MmcifConnectionError,MmcifConnectionErrorKind,ck::MmcifConnectionError,{
Cif=>ck::MmcifConnectionError::Cif(_),
SequenceId=>ck::MmcifConnectionError::SequenceId(_),
LabelAddress=>ck::MmcifConnectionError::LabelAddress(_),
MissingAddressIdentifiers=>ck::MmcifConnectionError::MissingAddressIdentifiers{..},
NoStructuralModels=>ck::MmcifConnectionError::NoStructuralModels,
ChainNameNotRepresentable=>ck::MmcifConnectionError::ChainNameNotRepresentable{..},
ResidueNameNotRepresentable=>ck::MmcifConnectionError::ResidueNameNotRepresentable{..},
ResidueAddressNotRepresentable=>ck::MmcifConnectionError::ResidueAddressNotRepresentable{..},
SymmetryOperatorOutsideSourceDefinedRange=>ck::MmcifConnectionError::SymmetryOperatorOutsideSourceDefinedRange{..},
SymmetryOperatorNotRepresentable=>ck::MmcifConnectionError::SymmetryOperatorNotRepresentable{..}
},{
#[wasm_bindgen(getter,js_name=partner)] pub fn partner(&self)->Option<usize>{match self.inner.as_ref(){
ck::MmcifConnectionError::MissingAddressIdentifiers{partner,..}=>Some(*partner),
ck::MmcifConnectionError::ChainNameNotRepresentable{partner,..}=>Some(*partner),
ck::MmcifConnectionError::ResidueNameNotRepresentable{partner,..}=>Some(*partner),
ck::MmcifConnectionError::ResidueAddressNotRepresentable{partner,..}=>Some(*partner),
_=>None,}}
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<MmcifConnectionErrorValueValue>{match self.inner.as_ref(){
ck::MmcifConnectionError::ChainNameNotRepresentable{value,..}=>Some(JsValue::from(value.to_string()).unchecked_into()),
ck::MmcifConnectionError::ResidueNameNotRepresentable{value,..}=>Some(JsValue::from(value.to_string()).unchecked_into()),
ck::MmcifConnectionError::ResidueAddressNotRepresentable{value,..}=>Some(JsValue::from(value.to_string()).unchecked_into()),
ck::MmcifConnectionError::SymmetryOperatorOutsideSourceDefinedRange{value,..}=>Some(JsValue::from(value.to_string()).unchecked_into()),
ck::MmcifConnectionError::SymmetryOperatorNotRepresentable{value,..}=>Some(JsValue::from(*value).unchecked_into()),
_=>None,}}
},source=MmcifConnectionErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifConnectionError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "string | number")]
    pub type MmcifConnectionErrorValueValue;
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | MmcifSequenceIdError")]
    #[derive(Clone)]
    pub type MmcifHelixErrorCause;
}
enum_error!(MmcifHelixError,MmcifHelixErrorKind,ck::MmcifHelixError,{
Cif=>ck::MmcifHelixError::Cif(_),
ChainNameNotRepresentable=>ck::MmcifHelixError::ChainNameNotRepresentable{..},
ResidueNameNotRepresentable=>ck::MmcifHelixError::ResidueNameNotRepresentable{..},
ResidueAddressNotRepresentable=>ck::MmcifHelixError::ResidueAddressNotRepresentable{..},
SequenceId=>ck::MmcifHelixError::SequenceId{..},
HelixClass=>ck::MmcifHelixError::HelixClass(_),
HelixLength=>ck::MmcifHelixError::HelixLength(_)
},{
#[wasm_bindgen(getter,js_name=endpoint)] pub fn endpoint(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifHelixError::ChainNameNotRepresentable{endpoint,..}=>Some(endpoint.to_string()),
ck::MmcifHelixError::ResidueNameNotRepresentable{endpoint,..}=>Some(endpoint.to_string()),
ck::MmcifHelixError::ResidueAddressNotRepresentable{endpoint,..}=>Some(endpoint.to_string()),
ck::MmcifHelixError::SequenceId{endpoint,..}=>Some(endpoint.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifHelixError::ChainNameNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifHelixError::ResidueNameNotRepresentable{value,..}=>Some(value.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=error)] pub fn error(&self)->Option<MmcifHelixErrorCause>{match self.inner.as_ref(){
ck::MmcifHelixError::SequenceId{error,..}=>Some(MmcifSequenceIdError::from_owner(self.inner.child(|e|match e{ck::MmcifHelixError::SequenceId{error,..}=>error,_=>unreachable!("immutable BIO error variant")})).unchecked_into()),
_=>None,}}
},source=MmcifHelixErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifHelixError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | MmcifSequenceIdError")]
    #[derive(Clone)]
    pub type MmcifSheetErrorCause;
}
enum_error!(MmcifSheetError,MmcifSheetErrorKind,ck::MmcifSheetError,{
Cif=>ck::MmcifSheetError::Cif(_),
ChainNameNotRepresentable=>ck::MmcifSheetError::ChainNameNotRepresentable{..},
ResidueNameNotRepresentable=>ck::MmcifSheetError::ResidueNameNotRepresentable{..},
ResidueAddressNotRepresentable=>ck::MmcifSheetError::ResidueAddressNotRepresentable{..},
SequenceId=>ck::MmcifSheetError::SequenceId{..}
},{
#[wasm_bindgen(getter,js_name=address)] pub fn address(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSheetError::ChainNameNotRepresentable{address,..}=>Some(address.to_string()),
ck::MmcifSheetError::ResidueNameNotRepresentable{address,..}=>Some(address.to_string()),
ck::MmcifSheetError::ResidueAddressNotRepresentable{address,..}=>Some(address.to_string()),
ck::MmcifSheetError::SequenceId{address,..}=>Some(address.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifSheetError::ChainNameNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifSheetError::ResidueNameNotRepresentable{value,..}=>Some(value.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=error)] pub fn error(&self)->Option<MmcifSheetErrorCause>{match self.inner.as_ref(){
ck::MmcifSheetError::SequenceId{error,..}=>Some(MmcifSequenceIdError::from_owner(self.inner.child(|e|match e{ck::MmcifSheetError::SequenceId{error,..}=>error,_=>unreachable!("immutable BIO error variant")})).unchecked_into()),
_=>None,}}
},source=MmcifSheetErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifSheetError>|native_source(owner).map(JsCast::unchecked_into));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | MmcifSequenceIdError")]
    #[derive(Clone)]
    pub type MmcifCisPepErrorCause;
}
enum_error!(MmcifCisPepError,MmcifCisPepErrorKind,ck::MmcifCisPepError,{
Cif=>ck::MmcifCisPepError::Cif(_),
ChainNameNotRepresentable=>ck::MmcifCisPepError::ChainNameNotRepresentable{..},
ResidueNameNotRepresentable=>ck::MmcifCisPepError::ResidueNameNotRepresentable{..},
ResidueAddressNotRepresentable=>ck::MmcifCisPepError::ResidueAddressNotRepresentable{..},
SequenceId=>ck::MmcifCisPepError::SequenceId{..}
},{
#[wasm_bindgen(getter,js_name=partner)] pub fn partner(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifCisPepError::ChainNameNotRepresentable{partner,..}=>Some(partner.to_string()),
ck::MmcifCisPepError::ResidueNameNotRepresentable{partner,..}=>Some(partner.to_string()),
ck::MmcifCisPepError::ResidueAddressNotRepresentable{partner,..}=>Some(partner.to_string()),
ck::MmcifCisPepError::SequenceId{partner,..}=>Some(partner.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=value)] pub fn value(&self)->Option<String>{match self.inner.as_ref(){
ck::MmcifCisPepError::ChainNameNotRepresentable{value,..}=>Some(value.to_string()),
ck::MmcifCisPepError::ResidueNameNotRepresentable{value,..}=>Some(value.to_string()),
_=>None,}}
#[wasm_bindgen(getter,js_name=error)] pub fn error(&self)->Option<MmcifCisPepErrorCause>{match self.inner.as_ref(){
ck::MmcifCisPepError::SequenceId{error,..}=>Some(MmcifSequenceIdError::from_owner(self.inner.child(|e|match e{ck::MmcifCisPepError::SequenceId{error,..}=>error,_=>unreachable!("immutable BIO error variant")})).unchecked_into()),
_=>None,}}
},source=MmcifCisPepErrorCause,cause=|owner:&ReadonlyValue<ck::MmcifCisPepError>|native_source(owner).map(JsCast::unchecked_into));

crate::enums::enum_value!(CifReadErrorKind,"CifReadErrorKind" {Syntax="syntax",MissingValue="missing_value",DuplicateName="duplicate_name",InvalidLoop="invalid_loop",InvalidValue="invalid_value",OutOfRange="out_of_range"});
crate::enums::enum_value!(AtomSiteIdentityField,"AtomSiteIdentityField" {AsymId="asym_id",CompId="comp_id",AtomId="atom_id",SeqId="seq_id"});
crate::enums::enum_value!(BioMmcifReadStage,"BioMmcifReadStage" {CifDocument="cif_document",CoordinateBlock="coordinate_block",CrystalCell="crystal_cell",Refinement="refinement",Tls="tls",Experimental="experimental",Reflections="reflections",Software="software",Ncs="ncs",FractionalTransform="fractional_transform",Origx="origx",AnisotropicU="anisotropic_u",AtomSites="atom_sites",EntitySequence="entity_sequence",Helices="helices",Sheets="sheets",Connections="connections",CisPeptides="cis_peptides",ModifiedResidues="modified_residues",Assemblies="assemblies",SiftsUnp="sifts_unp",CcdRestoration="ccd_restoration",Materialization="materialization",StructureValidation="structure_validation"});
error_projection!(CifReadError, ck::CifReadError, {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> CifReadErrorKind {
        self.inner.kind().into()
    }
    #[wasm_bindgen(getter)]
    pub fn source(&self) -> String {
        self.inner.source().into()
    }
    #[wasm_bindgen(getter)]
    pub fn line(&self) -> usize {
        self.inner.line()
    }
    #[wasm_bindgen(getter)]
    pub fn column(&self) -> usize {
        self.inner.column()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.message().into()
    }
});
error_projection!(
    MissingAtomSiteIdentityField,
    ck::MissingAtomSiteIdentityField,
    {
        #[wasm_bindgen(getter)]
        pub fn field(&self) -> AtomSiteIdentityField {
            self.inner.field().into()
        }
    }
);
error_projection!(
    MmcifAtomSiteRowOutOfBounds,
    ck::MmcifAtomSiteRowOutOfBounds,
    {
        #[wasm_bindgen(getter,js_name=rowIndex)]
        pub fn row_index(&self) -> usize {
            self.inner.row_index()
        }
        #[wasm_bindgen(getter,js_name=parsedCount)]
        pub fn parsed_count(&self) -> usize {
            self.inner.parsed_count()
        }
    }
);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "PdbResidueKeyError | CoordinateBlockSelectionError | BioStructureError | MmcifCrystalInfoError | MmcifSequenceIdError | MmcifTlsInfoError | MmcifSiftsUnpError | MmcifEntityInfoError | MmcifLabelAddressError | MmcifOperationExprError | MmcifAssemblyError | MmcifAtomSiteScalarError | MmcifModResError | MmcifCcdRestoreError | MmcifAtomSiteGroupingError | MmcifConnectionError | MmcifHelixError | MmcifSheetError | MmcifCisPepError | CifReadError | MissingAtomSiteIdentityField | MmcifAtomSiteRowOutOfBounds | FromUtf8Error"
    )]
    #[derive(Clone)]
    pub type MmcifCause;
}
error_projection!(
    BioMmcifReadError,
    ck::BioMmcifReadError,
    {
        #[wasm_bindgen(getter)]
        pub fn stage(&self) -> BioMmcifReadStage {
            self.inner.stage().into()
        }
    },
    source = MmcifCause,
    cause = |owner: &ReadonlyValue<ck::BioMmcifReadError>| native_source(owner)
        .map(JsCast::unchecked_into)
);
fn native_source<T: Error + 'static>(owner: &ReadonlyValue<T>) -> Option<JsValue> {
    Error::source(owner.as_ref())?;
    Some(project(owner.child(|v| {
        Error::source(v).expect("immutable native BIO source")
    })))
}
fn project(owner: ReadonlyValue<dyn Error>) -> JsValue {
    if owner.is::<ck::PdbResidueKeyError>() {
        return PdbResidueKeyError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::PdbResidueKeyError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::CoordinateBlockSelectionError>() {
        return CoordinateBlockSelectionError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::CoordinateBlockSelectionError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::BioStructureError>() {
        return BioStructureError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::BioStructureError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifCrystalInfoError>() {
        return MmcifCrystalInfoError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifCrystalInfoError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifSequenceIdError>() {
        return MmcifSequenceIdError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifSequenceIdError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifTlsInfoError>() {
        return MmcifTlsInfoError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifTlsInfoError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifSiftsUnpError>() {
        return MmcifSiftsUnpError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifSiftsUnpError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifEntityInfoError>() {
        return MmcifEntityInfoError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifEntityInfoError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifLabelAddressError>() {
        return MmcifLabelAddressError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifLabelAddressError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifOperationExprError>() {
        return MmcifOperationExprError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifOperationExprError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifAssemblyError>() {
        return MmcifAssemblyError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifAssemblyError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifAtomSiteScalarError>() {
        return MmcifAtomSiteScalarError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifAtomSiteScalarError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifModResError>() {
        return MmcifModResError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifModResError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifCcdRestoreError>() {
        return MmcifCcdRestoreError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifCcdRestoreError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifAtomSiteGroupingError>() {
        return MmcifAtomSiteGroupingError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifAtomSiteGroupingError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifConnectionError>() {
        return MmcifConnectionError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifConnectionError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifHelixError>() {
        return MmcifHelixError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifHelixError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifSheetError>() {
        return MmcifSheetError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifSheetError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifCisPepError>() {
        return MmcifCisPepError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifCisPepError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::CifReadError>() {
        return CifReadError::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::CifReadError>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MissingAtomSiteIdentityField>() {
        return MissingAtomSiteIdentityField::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MissingAtomSiteIdentityField>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<ck::MmcifAtomSiteRowOutOfBounds>() {
        return MmcifAtomSiteRowOutOfBounds::from_owner(owner.child(|v| {
            v.downcast_ref::<ck::MmcifAtomSiteRowOutOfBounds>()
                .expect("selected native BIO cause")
        }));
    }
    if owner.is::<std::string::FromUtf8Error>() {
        return FromUtf8Error::from_owner(owner.child(|v| {
            v.downcast_ref::<std::string::FromUtf8Error>()
                .expect("selected UTF-8 cause")
        }));
    }
    unreachable!("native mmCIF reader returned an undeclared cause")
}
