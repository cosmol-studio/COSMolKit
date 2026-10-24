//! PDB reader payloads retain exact native variants, not formatted substitutes.
use crate::bio_read_errors::{BioStructureError, PdbResidueKeyError};
use crate::errors::{enum_error, error_projection, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "BioStructureError")]
    pub type BioStructureErrorRef;
    #[wasm_bindgen(typescript_type = "PdbResidueKeyError")]
    pub type PdbResidueKeyErrorRef;
    #[wasm_bindgen(typescript_type = "number | null")]
    pub type NullableByte;
    #[wasm_bindgen(typescript_type = "PdbReaderCause")]
    #[derive(Clone)]
    pub type PdbReaderCauseRef;
}
crate::enums::enum_value!(PdbLinkNameCounter,"PdbLinkNameCounter"{Covale="covale",MetalC="metal_c"});
#[wasm_bindgen]
pub struct PdbDbRefError(ReadonlyValue<ck::PdbDbRefError>);
impl PdbDbRefError {
    fn from_owner(value: ReadonlyValue<ck::PdbDbRefError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbDbRefErrorKind {
    SourceIntegerOutsideDefinedRange,
    MissingCurrentDbRef,
}
#[wasm_bindgen]
impl PdbDbRefError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbDbRefErrorKind {
        match self.0.as_ref() {
            ck::PdbDbRefError::SourceIntegerOutsideDefinedRange { .. } => {
                PdbDbRefErrorKind::SourceIntegerOutsideDefinedRange
            }
            ck::PdbDbRefError::MissingCurrentDbRef => PdbDbRefErrorKind::MissingCurrentDbRef,
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbCrystalError(ReadonlyValue<ck::PdbCrystalError>);
impl PdbCrystalError {
    fn from_owner(value: ReadonlyValue<ck::PdbCrystalError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbCrystalErrorKind {
    BioStructureData,
    TextFieldNotUtf8,
}
#[wasm_bindgen]
impl PdbCrystalError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbCrystalErrorKind {
        match self.0.as_ref() {
            ck::PdbCrystalError::BioStructureData(_) => PdbCrystalErrorKind::BioStructureData,
            ck::PdbCrystalError::TextFieldNotUtf8 { .. } => PdbCrystalErrorKind::TextFieldNotUtf8,
        }
    }
    #[wasm_bindgen(getter,js_name=bioStructureData)]
    pub fn bio_structure_data(&self) -> Option<BioStructureErrorRef> {
        match self.0.as_ref() {
            ck::PdbCrystalError::BioStructureData(value) => Some(
                BioStructureError::from_owner(self.0.child(|e| match e {
                    ck::PdbCrystalError::BioStructureData(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbCrystalError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbCrystalError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbMtrixError(ReadonlyValue<ck::PdbMtrixError>);
impl PdbMtrixError {
    fn from_owner(value: ReadonlyValue<ck::PdbMtrixError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbMtrixErrorKind {
    TextFieldNotUtf8,
}
#[wasm_bindgen]
impl PdbMtrixError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbMtrixErrorKind {
        match self.0.as_ref() {
            ck::PdbMtrixError::TextFieldNotUtf8 { .. } => PdbMtrixErrorKind::TextFieldNotUtf8,
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbMtrixError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbMtrixError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbConectError(ReadonlyValue<ck::PdbConectError>);
impl PdbConectError {
    fn from_owner(value: ReadonlyValue<ck::PdbConectError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbConectErrorKind {
    SourceSerialOutsideDefinedRange,
}
#[wasm_bindgen]
impl PdbConectError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbConectErrorKind {
        match self.0.as_ref() {
            ck::PdbConectError::SourceSerialOutsideDefinedRange { .. } => {
                PdbConectErrorKind::SourceSerialOutsideDefinedRange
            }
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbConectError::SourceSerialOutsideDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbModResError(ReadonlyValue<ck::PdbModResError>);
impl PdbModResError {
    fn from_owner(value: ReadonlyValue<ck::PdbModResError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbModResErrorKind {
    ChainNameNotRepresentable,
    ResidueAddress,
    TextFieldNotUtf8,
}
#[wasm_bindgen]
impl PdbModResError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbModResErrorKind {
        match self.0.as_ref() {
            ck::PdbModResError::ChainNameNotRepresentable { .. } => {
                PdbModResErrorKind::ChainNameNotRepresentable
            }
            ck::PdbModResError::ResidueAddress(_) => PdbModResErrorKind::ResidueAddress,
            ck::PdbModResError::TextFieldNotUtf8 { .. } => PdbModResErrorKind::TextFieldNotUtf8,
        }
    }
    #[wasm_bindgen(getter,js_name=width)]
    pub fn width(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbModResError::ChainNameNotRepresentable { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=nonAscii)]
    pub fn non_ascii(&self) -> Option<NullableByte> {
        match self.0.as_ref() {
            ck::PdbModResError::ChainNameNotRepresentable { non_ascii, .. } => Some(
                non_ascii
                    .map_or(JsValue::NULL, JsValue::from)
                    .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=residueAddress)]
    pub fn residue_address(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbModResError::ResidueAddress(value) => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbModResError::ResidueAddress(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbModResError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbModResError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbHelixError(ReadonlyValue<ck::PdbHelixError>);
impl PdbHelixError {
    fn from_owner(value: ReadonlyValue<ck::PdbHelixError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbHelixErrorKind {
    ChainNameNotRepresentable,
    ResidueAddress,
    IntegerOutsideSourceDefinedRange,
}
#[wasm_bindgen]
impl PdbHelixError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbHelixErrorKind {
        match self.0.as_ref() {
            ck::PdbHelixError::ChainNameNotRepresentable { .. } => {
                PdbHelixErrorKind::ChainNameNotRepresentable
            }
            ck::PdbHelixError::ResidueAddress(_) => PdbHelixErrorKind::ResidueAddress,
            ck::PdbHelixError::IntegerOutsideSourceDefinedRange { .. } => {
                PdbHelixErrorKind::IntegerOutsideSourceDefinedRange
            }
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbHelixError::ChainNameNotRepresentable { field_offset, .. } => {
                Some(*field_offset)
            }
            ck::PdbHelixError::IntegerOutsideSourceDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=width)]
    pub fn width(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbHelixError::ChainNameNotRepresentable { width, .. } => Some(*width),
            ck::PdbHelixError::IntegerOutsideSourceDefinedRange { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=nonAscii)]
    pub fn non_ascii(&self) -> Option<NullableByte> {
        match self.0.as_ref() {
            ck::PdbHelixError::ChainNameNotRepresentable { non_ascii, .. } => Some(
                non_ascii
                    .map_or(JsValue::NULL, JsValue::from)
                    .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=residueAddress)]
    pub fn residue_address(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbHelixError::ResidueAddress(value) => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbHelixError::ResidueAddress(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbSheetError(ReadonlyValue<ck::PdbSheetError>);
impl PdbSheetError {
    fn from_owner(value: ReadonlyValue<ck::PdbSheetError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbSheetErrorKind {
    TextFieldNotUtf8,
    ChainNameNotRepresentable,
    ResidueAddress,
    IntegerOutsideSourceDefinedRange,
}
#[wasm_bindgen]
impl PdbSheetError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbSheetErrorKind {
        match self.0.as_ref() {
            ck::PdbSheetError::TextFieldNotUtf8 { .. } => PdbSheetErrorKind::TextFieldNotUtf8,
            ck::PdbSheetError::ChainNameNotRepresentable { .. } => {
                PdbSheetErrorKind::ChainNameNotRepresentable
            }
            ck::PdbSheetError::ResidueAddress { .. } => PdbSheetErrorKind::ResidueAddress,
            ck::PdbSheetError::IntegerOutsideSourceDefinedRange { .. } => {
                PdbSheetErrorKind::IntegerOutsideSourceDefinedRange
            }
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbSheetError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            ck::PdbSheetError::ChainNameNotRepresentable { field_offset, .. } => {
                Some(*field_offset)
            }
            ck::PdbSheetError::IntegerOutsideSourceDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbSheetError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=width)]
    pub fn width(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbSheetError::ChainNameNotRepresentable { width, .. } => Some(*width),
            ck::PdbSheetError::IntegerOutsideSourceDefinedRange { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=nonAscii)]
    pub fn non_ascii(&self) -> Option<NullableByte> {
        match self.0.as_ref() {
            ck::PdbSheetError::ChainNameNotRepresentable { non_ascii, .. } => Some(
                non_ascii
                    .map_or(JsValue::NULL, JsValue::from)
                    .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=sequenceOffset)]
    pub fn sequence_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbSheetError::ResidueAddress {
                sequence_offset, ..
            } => Some(*sequence_offset),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=source)]
    pub fn source(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbSheetError::ResidueAddress { source, .. } => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbSheetError::ResidueAddress { source, .. } => source,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbHeaderError(ReadonlyValue<ck::PdbHeaderError>);
impl PdbHeaderError {
    fn from_owner(value: ReadonlyValue<ck::PdbHeaderError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbHeaderErrorKind {
    TextFieldNotUtf8,
}
#[wasm_bindgen]
impl PdbHeaderError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbHeaderErrorKind {
        match self.0.as_ref() {
            ck::PdbHeaderError::TextFieldNotUtf8 { .. } => PdbHeaderErrorKind::TextFieldNotUtf8,
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbHeaderError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbHeaderError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbRemarkError(ReadonlyValue<ck::PdbRemarkError>);
impl PdbRemarkError {
    fn from_owner(value: ReadonlyValue<ck::PdbRemarkError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbRemarkErrorKind {
    LineLengthOutsideBuffer,
    TextFieldNotUtf8,
}
#[wasm_bindgen]
impl PdbRemarkError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbRemarkErrorKind {
        match self.0.as_ref() {
            ck::PdbRemarkError::LineLengthOutsideBuffer { .. } => {
                PdbRemarkErrorKind::LineLengthOutsideBuffer
            }
            ck::PdbRemarkError::TextFieldNotUtf8 { .. } => PdbRemarkErrorKind::TextFieldNotUtf8,
        }
    }
    #[wasm_bindgen(getter,js_name=lineLen)]
    pub fn line_len(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemarkError::LineLengthOutsideBuffer { line_len, .. } => Some(*line_len),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemarkError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbRemark3Error(ReadonlyValue<ck::PdbRemark3Error>);
impl PdbRemark3Error {
    fn from_owner(value: ReadonlyValue<ck::PdbRemark3Error>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbRemark3ErrorKind {
    IntegerOutsideSourceDefinedRange,
    TextNotUtf8,
    MissingContinuationTarget,
    TlsChainIdNotRepresentable,
}
#[wasm_bindgen]
impl PdbRemark3Error {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbRemark3ErrorKind {
        match self.0.as_ref() {
            ck::PdbRemark3Error::IntegerOutsideSourceDefinedRange => {
                PdbRemark3ErrorKind::IntegerOutsideSourceDefinedRange
            }
            ck::PdbRemark3Error::TextNotUtf8 { .. } => PdbRemark3ErrorKind::TextNotUtf8,
            ck::PdbRemark3Error::MissingContinuationTarget => {
                PdbRemark3ErrorKind::MissingContinuationTarget
            }
            ck::PdbRemark3Error::TlsChainIdNotRepresentable { .. } => {
                PdbRemark3ErrorKind::TlsChainIdNotRepresentable
            }
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark3Error::TextNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=width)]
    pub fn width(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark3Error::TlsChainIdNotRepresentable { width, .. } => Some(*width),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbRemark200Error(ReadonlyValue<ck::PdbRemark200Error>);
impl PdbRemark200Error {
    fn from_owner(value: ReadonlyValue<ck::PdbRemark200Error>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbRemark200ErrorKind {
    IntegerOutsideSourceDefinedRange,
    TextNotUtf8,
    MissingExperimentalRows,
    UndefinedPhNumericScreening,
    Software,
}
#[wasm_bindgen]
impl PdbRemark200Error {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbRemark200ErrorKind {
        match self.0.as_ref() {
            ck::PdbRemark200Error::IntegerOutsideSourceDefinedRange => {
                PdbRemark200ErrorKind::IntegerOutsideSourceDefinedRange
            }
            ck::PdbRemark200Error::TextNotUtf8 { .. } => PdbRemark200ErrorKind::TextNotUtf8,
            ck::PdbRemark200Error::MissingExperimentalRows => {
                PdbRemark200ErrorKind::MissingExperimentalRows
            }
            ck::PdbRemark200Error::UndefinedPhNumericScreening => {
                PdbRemark200ErrorKind::UndefinedPhNumericScreening
            }
            ck::PdbRemark200Error::Software(_) => PdbRemark200ErrorKind::Software,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark200Error::TextNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=software)]
    pub fn software(&self) -> Option<PdbRemark3Error> {
        match self.0.as_ref() {
            ck::PdbRemark200Error::Software(value) => {
                Some(PdbRemark3Error::from_owner(self.0.child(|e| match e {
                    ck::PdbRemark200Error::Software(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })))
            }
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbRemarkMetadataError(ReadonlyValue<ck::PdbRemarkMetadataError>);
impl PdbRemarkMetadataError {
    fn from_owner(value: ReadonlyValue<ck::PdbRemarkMetadataError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbRemarkMetadataErrorKind {
    IntegerOutsideSourceDefinedRange,
    UndefinedResolutionField,
    TextBoundaryNotUtf8,
}
#[wasm_bindgen]
impl PdbRemarkMetadataError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbRemarkMetadataErrorKind {
        match self.0.as_ref() {
            ck::PdbRemarkMetadataError::IntegerOutsideSourceDefinedRange => {
                PdbRemarkMetadataErrorKind::IntegerOutsideSourceDefinedRange
            }
            ck::PdbRemarkMetadataError::UndefinedResolutionField { .. } => {
                PdbRemarkMetadataErrorKind::UndefinedResolutionField
            }
            ck::PdbRemarkMetadataError::TextBoundaryNotUtf8 { .. } => {
                PdbRemarkMetadataErrorKind::TextBoundaryNotUtf8
            }
        }
    }
    #[wasm_bindgen(getter,js_name=remarkLength)]
    pub fn remark_length(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemarkMetadataError::UndefinedResolutionField { remark_length, .. } => {
                Some(*remark_length)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=offset)]
    pub fn offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemarkMetadataError::TextBoundaryNotUtf8 { offset, .. } => Some(*offset),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbRemark350Error(ReadonlyValue<ck::PdbRemark350Error>);
impl PdbRemark350Error {
    fn from_owner(value: ReadonlyValue<ck::PdbRemark350Error>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbRemark350ErrorKind {
    IntegerOutsideSourceDefinedRange,
    SourceUndefinedAccess,
    TextNotUtf8,
    StateInvariantLost,
}
#[wasm_bindgen]
impl PdbRemark350Error {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbRemark350ErrorKind {
        match self.0.as_ref() {
            ck::PdbRemark350Error::IntegerOutsideSourceDefinedRange => {
                PdbRemark350ErrorKind::IntegerOutsideSourceDefinedRange
            }
            ck::PdbRemark350Error::SourceUndefinedAccess { .. } => {
                PdbRemark350ErrorKind::SourceUndefinedAccess
            }
            ck::PdbRemark350Error::TextNotUtf8 { .. } => PdbRemark350ErrorKind::TextNotUtf8,
            ck::PdbRemark350Error::StateInvariantLost => PdbRemark350ErrorKind::StateInvariantLost,
        }
    }
    #[wasm_bindgen(getter,js_name=offset)]
    pub fn offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark350Error::SourceUndefinedAccess { offset, .. } => Some(*offset),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=width)]
    pub fn width(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark350Error::SourceUndefinedAccess { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=recordLength)]
    pub fn record_length(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark350Error::SourceUndefinedAccess { record_length, .. } => {
                Some(*record_length)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=fieldOffset)]
    pub fn field_offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark350Error::TextNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbRemark350Error::TextNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbLinkRecordError(ReadonlyValue<ck::PdbLinkRecordError>);
impl PdbLinkRecordError {
    fn from_owner(value: ReadonlyValue<ck::PdbLinkRecordError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbLinkRecordErrorKind {
    ChainNameNotAscii,
    AtomNameNotAscii,
    LinkIdNotUtf8,
    ResidueAddress,
    NameCounterOutsideSourceDefinedRange,
}
#[wasm_bindgen]
impl PdbLinkRecordError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbLinkRecordErrorKind {
        match self.0.as_ref() {
            ck::PdbLinkRecordError::ChainNameNotAscii { .. } => {
                PdbLinkRecordErrorKind::ChainNameNotAscii
            }
            ck::PdbLinkRecordError::AtomNameNotAscii { .. } => {
                PdbLinkRecordErrorKind::AtomNameNotAscii
            }
            ck::PdbLinkRecordError::LinkIdNotUtf8 { .. } => PdbLinkRecordErrorKind::LinkIdNotUtf8,
            ck::PdbLinkRecordError::ResidueAddress { .. } => PdbLinkRecordErrorKind::ResidueAddress,
            ck::PdbLinkRecordError::NameCounterOutsideSourceDefinedRange { .. } => {
                PdbLinkRecordErrorKind::NameCounterOutsideSourceDefinedRange
            }
        }
    }
    #[wasm_bindgen(getter,js_name=partner)]
    pub fn partner(&self) -> Option<u8> {
        match self.0.as_ref() {
            ck::PdbLinkRecordError::ChainNameNotAscii { partner, .. } => Some(*partner),
            ck::PdbLinkRecordError::AtomNameNotAscii { partner, .. } => Some(*partner),
            ck::PdbLinkRecordError::ResidueAddress { partner, .. } => Some(*partner),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=byte)]
    pub fn byte(&self) -> Option<u8> {
        match self.0.as_ref() {
            ck::PdbLinkRecordError::ChainNameNotAscii { byte, .. } => Some(*byte),
            ck::PdbLinkRecordError::AtomNameNotAscii { byte, .. } => Some(*byte),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbLinkRecordError::LinkIdNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=source)]
    pub fn source(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbLinkRecordError::ResidueAddress { source, .. } => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbLinkRecordError::ResidueAddress { source, .. } => source,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=counter)]
    pub fn counter(&self) -> Option<PdbLinkNameCounter> {
        match self.0.as_ref() {
            ck::PdbLinkRecordError::NameCounterOutsideSourceDefinedRange { counter, .. } => {
                Some((*counter).into())
            }
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbCisPepRecordError(ReadonlyValue<ck::PdbCisPepRecordError>);
impl PdbCisPepRecordError {
    fn from_owner(value: ReadonlyValue<ck::PdbCisPepRecordError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbCisPepRecordErrorKind {
    RequiredFieldsPastRecord,
    NoModelsAfterSourceDefault,
    SingleModelNumberUnavailable,
    ChainNameNotAscii,
    ResidueAddress,
    ModelNumberOutsideSourceDefinedRange,
}
#[wasm_bindgen]
impl PdbCisPepRecordError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbCisPepRecordErrorKind {
        match self.0.as_ref() {
            ck::PdbCisPepRecordError::RequiredFieldsPastRecord { .. } => {
                PdbCisPepRecordErrorKind::RequiredFieldsPastRecord
            }
            ck::PdbCisPepRecordError::NoModelsAfterSourceDefault => {
                PdbCisPepRecordErrorKind::NoModelsAfterSourceDefault
            }
            ck::PdbCisPepRecordError::SingleModelNumberUnavailable => {
                PdbCisPepRecordErrorKind::SingleModelNumberUnavailable
            }
            ck::PdbCisPepRecordError::ChainNameNotAscii { .. } => {
                PdbCisPepRecordErrorKind::ChainNameNotAscii
            }
            ck::PdbCisPepRecordError::ResidueAddress { .. } => {
                PdbCisPepRecordErrorKind::ResidueAddress
            }
            ck::PdbCisPepRecordError::ModelNumberOutsideSourceDefinedRange => {
                PdbCisPepRecordErrorKind::ModelNumberOutsideSourceDefinedRange
            }
        }
    }
    #[wasm_bindgen(getter,js_name=recordLen)]
    pub fn record_len(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbCisPepRecordError::RequiredFieldsPastRecord { record_len, .. } => {
                Some(*record_len)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=requiredLen)]
    pub fn required_len(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbCisPepRecordError::RequiredFieldsPastRecord { required_len, .. } => {
                Some(*required_len)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=partner)]
    pub fn partner(&self) -> Option<u8> {
        match self.0.as_ref() {
            ck::PdbCisPepRecordError::ChainNameNotAscii { partner, .. } => Some(*partner),
            ck::PdbCisPepRecordError::ResidueAddress { partner, .. } => Some(*partner),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=byte)]
    pub fn byte(&self) -> Option<u8> {
        match self.0.as_ref() {
            ck::PdbCisPepRecordError::ChainNameNotAscii { byte, .. } => Some(*byte),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=source)]
    pub fn source(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbCisPepRecordError::ResidueAddress { source, .. } => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbCisPepRecordError::ResidueAddress { source, .. } => source,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbModelTransitionError(ReadonlyValue<ck::PdbModelTransitionError>);
impl PdbModelTransitionError {
    fn from_owner(value: ReadonlyValue<ck::PdbModelTransitionError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbModelTransitionErrorKind {
    ModelWithoutEnd,
    DuplicateModelNumber,
    AtomBetweenModels,
    ModelNumberOutsideSourceIntRange,
    NonpolymerSubchainCounterOutsideSourceIntRange,
    AnisouWithoutAtom,
    DuplicateAnisou,
    AnisouIntegerOutsideSourceIntRange,
    ModelIdOutOfBounds,
    Structure,
}
#[wasm_bindgen]
impl PdbModelTransitionError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbModelTransitionErrorKind {
        match self.0.as_ref() {
            ck::PdbModelTransitionError::ModelWithoutEnd => {
                PdbModelTransitionErrorKind::ModelWithoutEnd
            }
            ck::PdbModelTransitionError::DuplicateModelNumber { .. } => {
                PdbModelTransitionErrorKind::DuplicateModelNumber
            }
            ck::PdbModelTransitionError::AtomBetweenModels => {
                PdbModelTransitionErrorKind::AtomBetweenModels
            }
            ck::PdbModelTransitionError::ModelNumberOutsideSourceIntRange => {
                PdbModelTransitionErrorKind::ModelNumberOutsideSourceIntRange
            }
            ck::PdbModelTransitionError::NonpolymerSubchainCounterOutsideSourceIntRange {
                ..
            } => PdbModelTransitionErrorKind::NonpolymerSubchainCounterOutsideSourceIntRange,
            ck::PdbModelTransitionError::AnisouWithoutAtom => {
                PdbModelTransitionErrorKind::AnisouWithoutAtom
            }
            ck::PdbModelTransitionError::DuplicateAnisou => {
                PdbModelTransitionErrorKind::DuplicateAnisou
            }
            ck::PdbModelTransitionError::AnisouIntegerOutsideSourceIntRange { .. } => {
                PdbModelTransitionErrorKind::AnisouIntegerOutsideSourceIntRange
            }
            ck::PdbModelTransitionError::ModelIdOutOfBounds { .. } => {
                PdbModelTransitionErrorKind::ModelIdOutOfBounds
            }
            ck::PdbModelTransitionError::Structure(_) => PdbModelTransitionErrorKind::Structure,
        }
    }
    #[wasm_bindgen(getter,js_name=number)]
    pub fn number(&self) -> Option<i32> {
        match self.0.as_ref() {
            ck::PdbModelTransitionError::DuplicateModelNumber { number, .. } => Some(*number),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=chain)]
    pub fn chain(&self) -> Option<String> {
        match self.0.as_ref() {
            ck::PdbModelTransitionError::NonpolymerSubchainCounterOutsideSourceIntRange {
                chain,
                ..
            } => Some(chain.as_str().into()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=offset)]
    pub fn offset(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbModelTransitionError::AnisouIntegerOutsideSourceIntRange { offset, .. } => {
                Some(*offset)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=id)]
    pub fn id(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbModelTransitionError::ModelIdOutOfBounds { id, .. } => Some(id.index()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=modelCount)]
    pub fn model_count(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbModelTransitionError::ModelIdOutOfBounds { model_count, .. } => {
                Some(*model_count)
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=structure)]
    pub fn structure(&self) -> Option<BioStructureErrorRef> {
        match self.0.as_ref() {
            ck::PdbModelTransitionError::Structure(value) => Some(
                BioStructureError::from_owner(self.0.child(|e| match e {
                    ck::PdbModelTransitionError::Structure(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbRemarkFinalizationError(ReadonlyValue<ck::PdbRemarkFinalizationError>);
impl PdbRemarkFinalizationError {
    fn from_owner(value: ReadonlyValue<ck::PdbRemarkFinalizationError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbRemarkFinalizationErrorKind {
    IntegerOutsideSourceDefinedRange,
    Remark3,
    Remark200,
    RemarkMetadata,
    Remark350,
}
#[wasm_bindgen]
impl PdbRemarkFinalizationError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbRemarkFinalizationErrorKind {
        match self.0.as_ref() {
            ck::PdbRemarkFinalizationError::IntegerOutsideSourceDefinedRange => {
                PdbRemarkFinalizationErrorKind::IntegerOutsideSourceDefinedRange
            }
            ck::PdbRemarkFinalizationError::Remark3(_) => PdbRemarkFinalizationErrorKind::Remark3,
            ck::PdbRemarkFinalizationError::Remark200(_) => {
                PdbRemarkFinalizationErrorKind::Remark200
            }
            ck::PdbRemarkFinalizationError::RemarkMetadata(_) => {
                PdbRemarkFinalizationErrorKind::RemarkMetadata
            }
            ck::PdbRemarkFinalizationError::Remark350(_) => {
                PdbRemarkFinalizationErrorKind::Remark350
            }
        }
    }
    #[wasm_bindgen(getter,js_name=remark3)]
    pub fn remark3(&self) -> Option<PdbRemark3Error> {
        match self.0.as_ref() {
            ck::PdbRemarkFinalizationError::Remark3(value) => {
                Some(PdbRemark3Error::from_owner(self.0.child(|e| match e {
                    ck::PdbRemarkFinalizationError::Remark3(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })))
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=remark200)]
    pub fn remark200(&self) -> Option<PdbRemark200Error> {
        match self.0.as_ref() {
            ck::PdbRemarkFinalizationError::Remark200(value) => {
                Some(PdbRemark200Error::from_owner(self.0.child(|e| match e {
                    ck::PdbRemarkFinalizationError::Remark200(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })))
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=remarkMetadata)]
    pub fn remark_metadata(&self) -> Option<PdbRemarkMetadataError> {
        match self.0.as_ref() {
            ck::PdbRemarkFinalizationError::RemarkMetadata(value) => Some(
                PdbRemarkMetadataError::from_owner(self.0.child(|e| match e {
                    ck::PdbRemarkFinalizationError::RemarkMetadata(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=remark350)]
    pub fn remark350(&self) -> Option<PdbRemark350Error> {
        match self.0.as_ref() {
            ck::PdbRemarkFinalizationError::Remark350(value) => {
                Some(PdbRemark350Error::from_owner(self.0.child(|e| match e {
                    ck::PdbRemarkFinalizationError::Remark350(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })))
            }
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbSsbondRecordError(ReadonlyValue<ck::PdbSsbondRecordError>);
impl PdbSsbondRecordError {
    fn from_owner(value: ReadonlyValue<ck::PdbSsbondRecordError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbSsbondRecordErrorKind {
    EmptyRecord,
    ChainNameNotAscii,
    ResidueAddress,
    CounterOutsideSourceIntRange,
}
#[wasm_bindgen]
impl PdbSsbondRecordError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbSsbondRecordErrorKind {
        match self.0.as_ref() {
            ck::PdbSsbondRecordError::EmptyRecord => PdbSsbondRecordErrorKind::EmptyRecord,
            ck::PdbSsbondRecordError::ChainNameNotAscii { .. } => {
                PdbSsbondRecordErrorKind::ChainNameNotAscii
            }
            ck::PdbSsbondRecordError::ResidueAddress { .. } => {
                PdbSsbondRecordErrorKind::ResidueAddress
            }
            ck::PdbSsbondRecordError::CounterOutsideSourceIntRange => {
                PdbSsbondRecordErrorKind::CounterOutsideSourceIntRange
            }
        }
    }
    #[wasm_bindgen(getter,js_name=partner)]
    pub fn partner(&self) -> Option<u8> {
        match self.0.as_ref() {
            ck::PdbSsbondRecordError::ChainNameNotAscii { partner, .. } => Some(*partner),
            ck::PdbSsbondRecordError::ResidueAddress { partner, .. } => Some(*partner),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=byte)]
    pub fn byte(&self) -> Option<u8> {
        match self.0.as_ref() {
            ck::PdbSsbondRecordError::ChainNameNotAscii { byte, .. } => Some(*byte),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=source)]
    pub fn source(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbSsbondRecordError::ResidueAddress { source, .. } => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbSsbondRecordError::ResidueAddress { source, .. } => source,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbConnectionFinalizationError(ReadonlyValue<ck::PdbConnectionFinalizationError>);
impl PdbConnectionFinalizationError {
    fn from_owner(value: ReadonlyValue<ck::PdbConnectionFinalizationError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbConnectionFinalizationErrorKind {
    SsbondRecord,
    LinkRecord,
    CisPepRecord,
    SsbondResolution,
}
#[wasm_bindgen]
impl PdbConnectionFinalizationError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbConnectionFinalizationErrorKind {
        match self.0.as_ref() {
            ck::PdbConnectionFinalizationError::SsbondRecord(_) => {
                PdbConnectionFinalizationErrorKind::SsbondRecord
            }
            ck::PdbConnectionFinalizationError::LinkRecord(_) => {
                PdbConnectionFinalizationErrorKind::LinkRecord
            }
            ck::PdbConnectionFinalizationError::CisPepRecord(_) => {
                PdbConnectionFinalizationErrorKind::CisPepRecord
            }
            ck::PdbConnectionFinalizationError::SsbondResolution(_) => {
                PdbConnectionFinalizationErrorKind::SsbondResolution
            }
        }
    }
    #[wasm_bindgen(getter,js_name=ssbondRecord)]
    pub fn ssbond_record(&self) -> Option<PdbSsbondRecordError> {
        match self.0.as_ref() {
            ck::PdbConnectionFinalizationError::SsbondRecord(value) => Some(
                PdbSsbondRecordError::from_owner(self.0.child(|e| match e {
                    ck::PdbConnectionFinalizationError::SsbondRecord(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=linkRecord)]
    pub fn link_record(&self) -> Option<PdbLinkRecordError> {
        match self.0.as_ref() {
            ck::PdbConnectionFinalizationError::LinkRecord(value) => {
                Some(PdbLinkRecordError::from_owner(self.0.child(|e| match e {
                    ck::PdbConnectionFinalizationError::LinkRecord(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })))
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=cisPepRecord)]
    pub fn cis_pep_record(&self) -> Option<PdbCisPepRecordError> {
        match self.0.as_ref() {
            ck::PdbConnectionFinalizationError::CisPepRecord(value) => Some(
                PdbCisPepRecordError::from_owner(self.0.child(|e| match e {
                    ck::PdbConnectionFinalizationError::CisPepRecord(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=ssbondResolution)]
    pub fn ssbond_resolution(&self) -> Option<PdbSsbondResolutionError> {
        match self.0.as_ref() {
            ck::PdbConnectionFinalizationError::SsbondResolution(value) => Some(
                PdbSsbondResolutionError::from_owner(self.0.child(|e| match e {
                    ck::PdbConnectionFinalizationError::SsbondResolution(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })),
            ),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbStructureFinalizationError(ReadonlyValue<ck::PdbStructureFinalizationError>);
impl PdbStructureFinalizationError {
    fn from_owner(value: ReadonlyValue<ck::PdbStructureFinalizationError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbStructureFinalizationErrorKind {
    Model,
    Connection,
    Remarks,
    Ccd,
    NoModelAfterDefault,
    SourceIdNotRepresentable,
    Utf8,
    Structure,
    StagingRowOrder,
}
#[wasm_bindgen]
impl PdbStructureFinalizationError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbStructureFinalizationErrorKind {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::Model(_) => PdbStructureFinalizationErrorKind::Model,
            ck::PdbStructureFinalizationError::Connection(_) => {
                PdbStructureFinalizationErrorKind::Connection
            }
            ck::PdbStructureFinalizationError::Remarks(_) => {
                PdbStructureFinalizationErrorKind::Remarks
            }
            ck::PdbStructureFinalizationError::Ccd(_) => PdbStructureFinalizationErrorKind::Ccd,
            ck::PdbStructureFinalizationError::NoModelAfterDefault => {
                PdbStructureFinalizationErrorKind::NoModelAfterDefault
            }
            ck::PdbStructureFinalizationError::SourceIdNotRepresentable { .. } => {
                PdbStructureFinalizationErrorKind::SourceIdNotRepresentable
            }
            ck::PdbStructureFinalizationError::Utf8 { .. } => {
                PdbStructureFinalizationErrorKind::Utf8
            }
            ck::PdbStructureFinalizationError::Structure(_) => {
                PdbStructureFinalizationErrorKind::Structure
            }
            ck::PdbStructureFinalizationError::StagingRowOrder { .. } => {
                PdbStructureFinalizationErrorKind::StagingRowOrder
            }
        }
    }
    #[wasm_bindgen(getter,js_name=model)]
    pub fn model(&self) -> Option<PdbModelTransitionError> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::Model(value) => Some(
                PdbModelTransitionError::from_owner(self.0.child(|e| match e {
                    ck::PdbStructureFinalizationError::Model(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=connection)]
    pub fn connection(&self) -> Option<PdbConnectionFinalizationError> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::Connection(value) => Some(
                PdbConnectionFinalizationError::from_owner(self.0.child(|e| match e {
                    ck::PdbStructureFinalizationError::Connection(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=remarks)]
    pub fn remarks(&self) -> Option<PdbRemarkFinalizationError> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::Remarks(value) => Some(
                PdbRemarkFinalizationError::from_owner(self.0.child(|e| match e {
                    ck::PdbStructureFinalizationError::Remarks(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=ccd)]
    pub fn ccd(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::Ccd(value) => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbStructureFinalizationError::Ccd(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=field)]
    pub fn field(&self) -> Option<String> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::SourceIdNotRepresentable { field, .. } => {
                Some(field.to_string())
            }
            ck::PdbStructureFinalizationError::Utf8 { field, .. } => Some(field.to_string()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::Utf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=structure)]
    pub fn structure(&self) -> Option<BioStructureErrorRef> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::Structure(value) => Some(
                BioStructureError::from_owner(self.0.child(|e| match e {
                    ck::PdbStructureFinalizationError::Structure(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=table)]
    pub fn table(&self) -> Option<String> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::StagingRowOrder { table, .. } => {
                Some(table.to_string())
            }
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=expected)]
    pub fn expected(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::StagingRowOrder { expected, .. } => Some(*expected),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=actual)]
    pub fn actual(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbStructureFinalizationError::StagingRowOrder { actual, .. } => Some(*actual),
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbAtomFieldError(ReadonlyValue<ck::PdbAtomFieldError>);
impl PdbAtomFieldError {
    fn from_owner(value: ReadonlyValue<ck::PdbAtomFieldError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbAtomFieldErrorKind {
    TooShort,
    SerialOutsideSourceDefinedRange,
    ChainNameNotRepresentable,
    AtomNameNotRepresentable,
    ResidueKey,
    Charge,
}
#[wasm_bindgen]
impl PdbAtomFieldError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbAtomFieldErrorKind {
        match self.0.as_ref() {
            ck::PdbAtomFieldError::TooShort { .. } => PdbAtomFieldErrorKind::TooShort,
            ck::PdbAtomFieldError::SerialOutsideSourceDefinedRange => {
                PdbAtomFieldErrorKind::SerialOutsideSourceDefinedRange
            }
            ck::PdbAtomFieldError::ChainNameNotRepresentable { .. } => {
                PdbAtomFieldErrorKind::ChainNameNotRepresentable
            }
            ck::PdbAtomFieldError::AtomNameNotRepresentable { .. } => {
                PdbAtomFieldErrorKind::AtomNameNotRepresentable
            }
            ck::PdbAtomFieldError::ResidueKey(_) => PdbAtomFieldErrorKind::ResidueKey,
            ck::PdbAtomFieldError::Charge(_) => PdbAtomFieldErrorKind::Charge,
        }
    }
    #[wasm_bindgen(getter,js_name=len)]
    pub fn len(&self) -> Option<usize> {
        match self.0.as_ref() {
            ck::PdbAtomFieldError::TooShort { len, .. } => Some(*len),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=bytes)]
    pub fn bytes(&self) -> Option<Vec<u8>> {
        match self.0.as_ref() {
            ck::PdbAtomFieldError::ChainNameNotRepresentable { bytes, .. } => Some(bytes.to_vec()),
            ck::PdbAtomFieldError::AtomNameNotRepresentable { bytes, .. } => Some(bytes.to_vec()),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=residueKey)]
    pub fn residue_key(&self) -> Option<PdbResidueKeyErrorRef> {
        match self.0.as_ref() {
            ck::PdbAtomFieldError::ResidueKey(value) => Some(
                PdbResidueKeyError::from_owner(self.0.child(|e| match e {
                    ck::PdbAtomFieldError::ResidueKey(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                }))
                .unchecked_into(),
            ),
            _ => None,
        }
    }
    #[wasm_bindgen(getter,js_name=charge)]
    pub fn charge(&self) -> Option<PdbChargeError> {
        match self.0.as_ref() {
            ck::PdbAtomFieldError::Charge(value) => {
                Some(PdbChargeError::from_owner(self.0.child(|e| match e {
                    ck::PdbAtomFieldError::Charge(value) => value,
                    _ => unreachable!("immutable PDB error payload"),
                })))
            }
            _ => None,
        }
    }
}
#[wasm_bindgen]
pub struct PdbInputFormatError(ReadonlyValue<ck::PdbInputFormatError>);
impl PdbInputFormatError {
    fn from_owner(value: ReadonlyValue<ck::PdbInputFormatError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
pub enum PdbInputFormatErrorKind {
    Cif,
    Mmjson,
}
#[wasm_bindgen]
impl PdbInputFormatError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> PdbInputFormatErrorKind {
        match self.0.as_ref() {
            ck::PdbInputFormatError::Cif { .. } => PdbInputFormatErrorKind::Cif,
            ck::PdbInputFormatError::Mmjson { .. } => PdbInputFormatErrorKind::Mmjson,
        }
    }
    #[wasm_bindgen(getter,js_name=source)]
    pub fn source(&self) -> Option<String> {
        match self.0.as_ref() {
            ck::PdbInputFormatError::Cif { source, .. } => Some(source.to_string()),
            ck::PdbInputFormatError::Mmjson { source, .. } => Some(source.to_string()),
            _ => None,
        }
    }
}
enum_error!(PdbReaderCause,PdbReaderCauseKind,ck::PdbReaderCause,{LineNumberOutsideSourceIntRange=>ck::PdbReaderCause::LineNumberOutsideSourceIntRange,InputFormat=>ck::PdbReaderCause::InputFormat(_),AtomFields=>ck::PdbReaderCause::AtomFields(_),ModelTransition=>ck::PdbReaderCause::ModelTransition(_),Remark=>ck::PdbReaderCause::Remark(_),Conect=>ck::PdbReaderCause::Conect(_),ModRes=>ck::PdbReaderCause::ModRes(_),Helix=>ck::PdbReaderCause::Helix(_),Sheet=>ck::PdbReaderCause::Sheet(_),Header=>ck::PdbReaderCause::Header(_),DbRef=>ck::PdbReaderCause::DbRef(_),Mtrix=>ck::PdbReaderCause::Mtrix(_),Crystal=>ck::PdbReaderCause::Crystal(_),Finalization=>ck::PdbReaderCause::Finalization(_)},{
#[wasm_bindgen(getter,js_name=inputFormat)] pub fn input_format(&self)->Option<PdbInputFormatError>{match self.inner.as_ref(){
ck::PdbReaderCause::InputFormat(value)=>Some(PdbInputFormatError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::InputFormat(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=atomFields)] pub fn atom_fields(&self)->Option<PdbAtomFieldError>{match self.inner.as_ref(){
ck::PdbReaderCause::AtomFields(value)=>Some(PdbAtomFieldError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::AtomFields(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=modelTransition)] pub fn model_transition(&self)->Option<PdbModelTransitionError>{match self.inner.as_ref(){
ck::PdbReaderCause::ModelTransition(value)=>Some(PdbModelTransitionError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::ModelTransition(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=remark)] pub fn remark(&self)->Option<PdbRemarkError>{match self.inner.as_ref(){
ck::PdbReaderCause::Remark(value)=>Some(PdbRemarkError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Remark(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=conect)] pub fn conect(&self)->Option<PdbConectError>{match self.inner.as_ref(){
ck::PdbReaderCause::Conect(value)=>Some(PdbConectError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Conect(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=modRes)] pub fn mod_res(&self)->Option<PdbModResError>{match self.inner.as_ref(){
ck::PdbReaderCause::ModRes(value)=>Some(PdbModResError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::ModRes(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=helix)] pub fn helix(&self)->Option<PdbHelixError>{match self.inner.as_ref(){
ck::PdbReaderCause::Helix(value)=>Some(PdbHelixError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Helix(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=sheet)] pub fn sheet(&self)->Option<PdbSheetError>{match self.inner.as_ref(){
ck::PdbReaderCause::Sheet(value)=>Some(PdbSheetError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Sheet(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=header)] pub fn header(&self)->Option<PdbHeaderError>{match self.inner.as_ref(){
ck::PdbReaderCause::Header(value)=>Some(PdbHeaderError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Header(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=dbRef)] pub fn db_ref(&self)->Option<PdbDbRefError>{match self.inner.as_ref(){
ck::PdbReaderCause::DbRef(value)=>Some(PdbDbRefError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::DbRef(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=mtrix)] pub fn mtrix(&self)->Option<PdbMtrixError>{match self.inner.as_ref(){
ck::PdbReaderCause::Mtrix(value)=>Some(PdbMtrixError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Mtrix(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=crystal)] pub fn crystal(&self)->Option<PdbCrystalError>{match self.inner.as_ref(){
ck::PdbReaderCause::Crystal(value)=>Some(PdbCrystalError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Crystal(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
#[wasm_bindgen(getter,js_name=finalization)] pub fn finalization(&self)->Option<PdbStructureFinalizationError>{match self.inner.as_ref(){
ck::PdbReaderCause::Finalization(value)=>Some(PdbStructureFinalizationError::from_owner(self.inner.child(|e|match e{ck::PdbReaderCause::Finalization(value)=>value,_=>unreachable!("immutable PDB error payload")}))),
_=>None,}}
});

#[wasm_bindgen]
pub struct PdbSsbondResolutionError(ReadonlyValue<ck::PdbSsbondResolutionError>);
impl PdbSsbondResolutionError {
    fn from_owner(value: ReadonlyValue<ck::PdbSsbondResolutionError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
impl PdbSsbondResolutionError {
    #[wasm_bindgen(getter)]
    pub fn partner(&self) -> u8 {
        self.0.partner()
    }
    #[wasm_bindgen(getter,js_name=atomName)]
    pub fn atom_name(&self) -> String {
        self.0.atom_name().into()
    }
}
#[wasm_bindgen]
pub struct PdbChargeError(ReadonlyValue<ck::PdbChargeError>);
impl PdbChargeError {
    fn from_owner(value: ReadonlyValue<ck::PdbChargeError>) -> Self {
        Self(value)
    }
}
#[wasm_bindgen]
impl PdbChargeError {
    #[wasm_bindgen(getter)]
    pub fn digit(&self) -> u8 {
        self.0.digit()
    }
    #[wasm_bindgen(getter)]
    pub fn sign(&self) -> u8 {
        self.0.sign()
    }
}
crate::enums::enum_value!(BioPdbReadStage,"BioPdbReadStage" {Stream="stream",Record="record",Finalization="finalization"});
error_projection!(
    BioPdbReadError,
    ck::BioPdbReadError,
    {
        #[wasm_bindgen(getter,js_name=lineNumber)]
        pub fn line_number(&self) -> Option<i32> {
            self.inner.line_number()
        }
        #[wasm_bindgen(getter)]
        pub fn stage(&self) -> BioPdbReadStage {
            self.inner.stage().into()
        }
        #[wasm_bindgen(getter,js_name=recordTag)]
        pub fn record_tag(&self) -> Option<Vec<u8>> {
            self.inner.record_tag().map(|v| v.to_vec())
        }
    },
    source = PdbReaderCauseRef,
    cause = |owner: &ReadonlyValue<ck::BioPdbReadError>| Some(
        PdbReaderCause::from_owner(owner.child(|v| {
            Error::source(v)
                .expect("native PDB cause")
                .downcast_ref::<ck::PdbReaderCause>()
                .expect("native PDB cause type")
        }))
        .unchecked_into()
    )
);
