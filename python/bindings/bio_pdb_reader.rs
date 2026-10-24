//! Concrete native PDB reader errors, readonly payload values and thin constructors.
use crate::bio_errors::BioStructureError;
use crate::bio_mmcif_reader::PdbResidueKeyError;
use crate::bio_structure::BioStructure;
use crate::chemical_values::enum_value;
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::ErrorRef;
use crate::public_selectors::BioPdbReadStage;
use crate::query_payloads::value_kind;
use crate::retained_error::RetainedError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;
enum_value!(PdbLinkNameCounter,{Covale,MetalC});

value_kind!(PdbDbRefErrorKind,PdbDbRefError,{SourceIntegerOutsideDefinedRange {field_offset},MissingCurrentDbRef});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbDbRefError {
    inner: ck::PdbDbRefError,
}
#[pymethods]
impl PdbDbRefError {
    #[getter]
    fn variant(&self) -> PdbDbRefErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbDbRefError::SourceIntegerOutsideDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
}
value_kind!(PdbMtrixErrorKind, PdbMtrixError, {
    TextFieldNotUtf8 {
        field_offset,
        valid_up_to,
    }
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbMtrixError {
    inner: ck::PdbMtrixError,
}
#[pymethods]
impl PdbMtrixError {
    #[getter]
    fn variant(&self) -> PdbMtrixErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbMtrixError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbMtrixError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
value_kind!(PdbConectErrorKind, PdbConectError, {
    SourceSerialOutsideDefinedRange { field_offset }
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbConectError {
    inner: ck::PdbConectError,
}
#[pymethods]
impl PdbConectError {
    #[getter]
    fn variant(&self) -> PdbConectErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbConectError::SourceSerialOutsideDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
}
value_kind!(PdbHelixErrorKind,PdbHelixError,{ChainNameNotRepresentable {field_offset,width,non_ascii},ResidueAddress(_),IntegerOutsideSourceDefinedRange {field_offset,width}});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbHelixError {
    inner: ck::PdbHelixError,
}
#[pymethods]
impl PdbHelixError {
    #[getter]
    fn variant(&self) -> PdbHelixErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbHelixError::ChainNameNotRepresentable { field_offset, .. }
            | ck::PdbHelixError::IntegerOutsideSourceDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
    #[getter]
    fn width(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbHelixError::ChainNameNotRepresentable { width, .. }
            | ck::PdbHelixError::IntegerOutsideSourceDefinedRange { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[getter]
    fn non_ascii(&self) -> Option<Option<u8>> {
        match &self.inner {
            ck::PdbHelixError::ChainNameNotRepresentable { non_ascii, .. } => Some(*non_ascii),
            _ => None,
        }
    }
    #[getter]
    fn residue_address(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match &self.inner {
            ck::PdbHelixError::ResidueAddress(source) => PdbResidueKeyError::from_rust(*source)
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
}
value_kind!(PdbModResErrorKind,PdbModResError,{ChainNameNotRepresentable {width,non_ascii},ResidueAddress(_),TextFieldNotUtf8 {field_offset,valid_up_to}});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbModResError {
    inner: ck::PdbModResError,
}
#[pymethods]
impl PdbModResError {
    #[getter]
    fn variant(&self) -> PdbModResErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn width(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbModResError::ChainNameNotRepresentable { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[getter]
    fn non_ascii(&self) -> Option<Option<u8>> {
        match &self.inner {
            ck::PdbModResError::ChainNameNotRepresentable { non_ascii, .. } => Some(*non_ascii),
            _ => None,
        }
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbModResError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbModResError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[getter]
    fn residue_address(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match &self.inner {
            ck::PdbModResError::ResidueAddress(source) => PdbResidueKeyError::from_rust(*source)
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
}
value_kind!(PdbSheetErrorKind,PdbSheetError,{TextFieldNotUtf8 {field_offset,valid_up_to},ChainNameNotRepresentable {field_offset,width,non_ascii},ResidueAddress {sequence_offset,source},IntegerOutsideSourceDefinedRange {field_offset,width}});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbSheetError {
    inner: ck::PdbSheetError,
}
#[pymethods]
impl PdbSheetError {
    #[getter]
    fn variant(&self) -> PdbSheetErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbSheetError::TextFieldNotUtf8 { field_offset, .. }
            | ck::PdbSheetError::ChainNameNotRepresentable { field_offset, .. }
            | ck::PdbSheetError::IntegerOutsideSourceDefinedRange { field_offset, .. } => {
                Some(*field_offset)
            }
            _ => None,
        }
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbSheetError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[getter]
    fn width(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbSheetError::ChainNameNotRepresentable { width, .. }
            | ck::PdbSheetError::IntegerOutsideSourceDefinedRange { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[getter]
    fn non_ascii(&self) -> Option<Option<u8>> {
        match &self.inner {
            ck::PdbSheetError::ChainNameNotRepresentable { non_ascii, .. } => Some(*non_ascii),
            _ => None,
        }
    }
    #[getter]
    fn sequence_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbSheetError::ResidueAddress {
                sequence_offset, ..
            } => Some(*sequence_offset),
            _ => None,
        }
    }
    #[getter]
    fn source(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match &self.inner {
            ck::PdbSheetError::ResidueAddress { source, .. } => {
                PdbResidueKeyError::from_rust(*source)
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
}
value_kind!(PdbRemark3ErrorKind,PdbRemark3Error,{IntegerOutsideSourceDefinedRange,TextNotUtf8 {valid_up_to},MissingContinuationTarget,TlsChainIdNotRepresentable {width}});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbRemark3Error {
    inner: ck::PdbRemark3Error,
}
#[pymethods]
impl PdbRemark3Error {
    #[getter]
    fn variant(&self) -> PdbRemark3ErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark3Error::TextNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[getter]
    fn width(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark3Error::TlsChainIdNotRepresentable { width, .. } => Some(*width),
            _ => None,
        }
    }
}
value_kind!(PdbRemark200ErrorKind,PdbRemark200Error,{IntegerOutsideSourceDefinedRange,TextNotUtf8 {valid_up_to},MissingExperimentalRows,UndefinedPhNumericScreening,Software(_)});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbRemark200Error {
    inner: ck::PdbRemark200Error,
}
#[pymethods]
impl PdbRemark200Error {
    #[getter]
    fn variant(&self) -> PdbRemark200ErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark200Error::TextNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[getter]
    fn software(&self) -> Option<PdbRemark3Error> {
        match self.inner {
            ck::PdbRemark200Error::Software(inner) => Some(PdbRemark3Error { inner }),
            _ => None,
        }
    }
}
value_kind!(PdbRemarkMetadataErrorKind,PdbRemarkMetadataError,{IntegerOutsideSourceDefinedRange,UndefinedResolutionField {remark_length},TextBoundaryNotUtf8 {offset}});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbRemarkMetadataError {
    inner: ck::PdbRemarkMetadataError,
}
#[pymethods]
impl PdbRemarkMetadataError {
    #[getter]
    fn variant(&self) -> PdbRemarkMetadataErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn remark_length(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemarkMetadataError::UndefinedResolutionField { remark_length, .. } => {
                Some(*remark_length)
            }
            _ => None,
        }
    }
    #[getter]
    fn offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemarkMetadataError::TextBoundaryNotUtf8 { offset, .. } => Some(*offset),
            _ => None,
        }
    }
}
value_kind!(PdbRemark350ErrorKind,PdbRemark350Error,{IntegerOutsideSourceDefinedRange,SourceUndefinedAccess {offset,width,record_length},TextNotUtf8 {field_offset,valid_up_to},StateInvariantLost});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbRemark350Error {
    inner: ck::PdbRemark350Error,
}
#[pymethods]
impl PdbRemark350Error {
    #[getter]
    fn variant(&self) -> PdbRemark350ErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark350Error::SourceUndefinedAccess { offset, .. } => Some(*offset),
            _ => None,
        }
    }
    #[getter]
    fn width(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark350Error::SourceUndefinedAccess { width, .. } => Some(*width),
            _ => None,
        }
    }
    #[getter]
    fn record_length(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark350Error::SourceUndefinedAccess { record_length, .. } => {
                Some(*record_length)
            }
            _ => None,
        }
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark350Error::TextNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbRemark350Error::TextNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
}
value_kind!(PdbSsbondRecordErrorKind,PdbSsbondRecordError,{EmptyRecord,ChainNameNotAscii {partner,byte},ResidueAddress {partner,source},CounterOutsideSourceIntRange});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbSsbondRecordError {
    inner: ck::PdbSsbondRecordError,
}
#[pymethods]
impl PdbSsbondRecordError {
    #[getter]
    fn variant(&self) -> PdbSsbondRecordErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn partner(&self) -> Option<u8> {
        match &self.inner {
            ck::PdbSsbondRecordError::ChainNameNotAscii { partner, .. }
            | ck::PdbSsbondRecordError::ResidueAddress { partner, .. } => Some(*partner),
            _ => None,
        }
    }
    #[getter]
    fn byte(&self) -> Option<u8> {
        match &self.inner {
            ck::PdbSsbondRecordError::ChainNameNotAscii { byte, .. } => Some(*byte),
            _ => None,
        }
    }
    #[getter]
    fn source(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match &self.inner {
            ck::PdbSsbondRecordError::ResidueAddress { source, .. } => {
                PdbResidueKeyError::from_rust(*source)
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
}
value_kind!(PdbLinkRecordErrorKind,PdbLinkRecordError,{ChainNameNotAscii {partner,byte},AtomNameNotAscii {partner,byte},LinkIdNotUtf8 {valid_up_to},ResidueAddress {partner,source},NameCounterOutsideSourceDefinedRange {counter}});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbLinkRecordError {
    inner: ck::PdbLinkRecordError,
}
#[pymethods]
impl PdbLinkRecordError {
    #[getter]
    fn variant(&self) -> PdbLinkRecordErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn partner(&self) -> Option<u8> {
        match &self.inner {
            ck::PdbLinkRecordError::ChainNameNotAscii { partner, .. }
            | ck::PdbLinkRecordError::AtomNameNotAscii { partner, .. }
            | ck::PdbLinkRecordError::ResidueAddress { partner, .. } => Some(*partner),
            _ => None,
        }
    }
    #[getter]
    fn byte(&self) -> Option<u8> {
        match &self.inner {
            ck::PdbLinkRecordError::ChainNameNotAscii { byte, .. }
            | ck::PdbLinkRecordError::AtomNameNotAscii { byte, .. } => Some(*byte),
            _ => None,
        }
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbLinkRecordError::LinkIdNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[getter]
    fn source(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match &self.inner {
            ck::PdbLinkRecordError::ResidueAddress { source, .. } => {
                PdbResidueKeyError::from_rust(*source)
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
    #[getter]
    fn counter(&self) -> Option<PdbLinkNameCounter> {
        match self.inner {
            ck::PdbLinkRecordError::NameCounterOutsideSourceDefinedRange { counter } => {
                Some(counter.into())
            }
            _ => None,
        }
    }
}
value_kind!(PdbCisPepRecordErrorKind,PdbCisPepRecordError,{RequiredFieldsPastRecord {record_len,required_len},NoModelsAfterSourceDefault,SingleModelNumberUnavailable,ChainNameNotAscii {partner,byte},ResidueAddress {partner,source},ModelNumberOutsideSourceDefinedRange});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbCisPepRecordError {
    inner: ck::PdbCisPepRecordError,
}
#[pymethods]
impl PdbCisPepRecordError {
    #[getter]
    fn variant(&self) -> PdbCisPepRecordErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn record_len(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbCisPepRecordError::RequiredFieldsPastRecord { record_len, .. } => {
                Some(*record_len)
            }
            _ => None,
        }
    }
    #[getter]
    fn required_len(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbCisPepRecordError::RequiredFieldsPastRecord { required_len, .. } => {
                Some(*required_len)
            }
            _ => None,
        }
    }
    #[getter]
    fn partner(&self) -> Option<u8> {
        match &self.inner {
            ck::PdbCisPepRecordError::ChainNameNotAscii { partner, .. }
            | ck::PdbCisPepRecordError::ResidueAddress { partner, .. } => Some(*partner),
            _ => None,
        }
    }
    #[getter]
    fn byte(&self) -> Option<u8> {
        match &self.inner {
            ck::PdbCisPepRecordError::ChainNameNotAscii { byte, .. } => Some(*byte),
            _ => None,
        }
    }
    #[getter]
    fn source(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match &self.inner {
            ck::PdbCisPepRecordError::ResidueAddress { source, .. } => {
                PdbResidueKeyError::from_rust(*source)
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
}
value_kind!(PdbCrystalErrorKind,PdbCrystalError,{BioStructureData(_),TextFieldNotUtf8 {field_offset,valid_up_to}});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbCrystalError {
    inner: ck::PdbCrystalError,
}
#[pymethods]
impl PdbCrystalError {
    #[getter]
    fn variant(&self) -> PdbCrystalErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn field_offset(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbCrystalError::TextFieldNotUtf8 { field_offset, .. } => Some(*field_offset),
            _ => None,
        }
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner {
            ck::PdbCrystalError::TextFieldNotUtf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[getter]
    fn bio_structure_data(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<BioStructureError>>> {
        match &self.inner {
            ck::PdbCrystalError::BioStructureData(inner) => {
                BioStructureError::from_rust(inner.clone())
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
}
value_kind!(PdbConnectionFinalizationErrorKind,PdbConnectionFinalizationError,{SsbondRecord(_),LinkRecord(_),CisPepRecord(_),SsbondResolution(_)});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbConnectionFinalizationError {
    inner: ck::PdbConnectionFinalizationError,
}
#[pymethods]
impl PdbConnectionFinalizationError {
    #[getter]
    fn variant(&self) -> PdbConnectionFinalizationErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn ssbond_record(&self) -> Option<PdbSsbondRecordError> {
        match &self.inner {
            ck::PdbConnectionFinalizationError::SsbondRecord(inner) => Some(PdbSsbondRecordError {
                inner: inner.clone(),
            }),
            _ => None,
        }
    }
    #[getter]
    fn link_record(&self) -> Option<PdbLinkRecordError> {
        match &self.inner {
            ck::PdbConnectionFinalizationError::LinkRecord(inner) => Some(PdbLinkRecordError {
                inner: inner.clone(),
            }),
            _ => None,
        }
    }
    #[getter]
    fn cis_pep_record(&self) -> Option<PdbCisPepRecordError> {
        match &self.inner {
            ck::PdbConnectionFinalizationError::CisPepRecord(inner) => Some(PdbCisPepRecordError {
                inner: inner.clone(),
            }),
            _ => None,
        }
    }
    #[getter]
    fn ssbond_resolution(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbSsbondResolutionError>>> {
        match &self.inner {
            ck::PdbConnectionFinalizationError::SsbondResolution(inner) => {
                PdbSsbondResolutionError::from_rust(inner.clone())
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
}
value_kind!(PdbRemarkFinalizationErrorKind,PdbRemarkFinalizationError,{IntegerOutsideSourceDefinedRange,Remark3(_),Remark200(_),RemarkMetadata(_),Remark350(_)});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbRemarkFinalizationError {
    inner: ck::PdbRemarkFinalizationError,
}
#[pymethods]
impl PdbRemarkFinalizationError {
    #[getter]
    fn variant(&self) -> PdbRemarkFinalizationErrorKind {
        (&self.inner).into()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
    #[getter]
    fn remark3(&self) -> Option<PdbRemark3Error> {
        match self.inner {
            ck::PdbRemarkFinalizationError::Remark3(inner) => Some(PdbRemark3Error { inner }),
            _ => None,
        }
    }
    #[getter]
    fn remark200(&self) -> Option<PdbRemark200Error> {
        match self.inner {
            ck::PdbRemarkFinalizationError::Remark200(inner) => Some(PdbRemark200Error { inner }),
            _ => None,
        }
    }
    #[getter]
    fn remark_metadata(&self) -> Option<PdbRemarkMetadataError> {
        match self.inner {
            ck::PdbRemarkFinalizationError::RemarkMetadata(inner) => {
                Some(PdbRemarkMetadataError { inner })
            }
            _ => None,
        }
    }
    #[getter]
    fn remark350(&self) -> Option<PdbRemark350Error> {
        match self.inner {
            ck::PdbRemarkFinalizationError::Remark350(inner) => Some(PdbRemark350Error { inner }),
            _ => None,
        }
    }
}

error_projection!(PdbRemarkError,PdbRemarkErrorKind,ck::PdbRemarkError,"bio",{
LineLengthOutsideBuffer{line_len:usize},TextFieldNotUtf8{valid_up_to:usize},
},cause=|_:&ck::PdbRemarkError|None::<PyErr>);
error_projection!(PdbHeaderError,PdbHeaderErrorKind,ck::PdbHeaderError,"bio",{
TextFieldNotUtf8{field_offset:usize,valid_up_to:usize},
},cause=|_:&ck::PdbHeaderError|None::<PyErr>);
error_projection!(PdbInputFormatError,PdbInputFormatErrorKind,ck::PdbInputFormatError,"bio",{
Cif{source},Mmjson{source},
},cause=|_:&ck::PdbInputFormatError|None::<PyErr>,methods={
#[getter]
fn source(&self)->&str{match &self.inner{ck::PdbInputFormatError::Cif{source}|ck::PdbInputFormatError::Mmjson{source}=>source}}
});
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct PdbSsbondResolutionError {
    inner: ck::PdbSsbondResolutionError,
}
}
impl PdbSsbondResolutionError {
    pub(crate) fn from_rust(inner: ck::PdbSsbondResolutionError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}
crate::error_framework::error_projection! {methods impl PdbSsbondResolutionError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn partner(&self) -> u8 {
        self.inner.partner()
    }
    #[getter]
    fn atom_name(&self) -> &str {
        self.inner.atom_name()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
error_projection!(retained PdbModelTransitionError,PdbModelTransitionErrorKind,ck::PdbModelTransitionError,"bio",{
ModelWithoutEnd,DuplicateModelNumber{number:i32},AtomBetweenModels,ModelNumberOutsideSourceIntRange,
NonpolymerSubchainCounterOutsideSourceIntRange{chain},AnisouWithoutAtom,DuplicateAnisou,
AnisouIntegerOutsideSourceIntRange{offset:usize},ModelIdOutOfBounds{id,model_count:usize},Structure(_),
},cause=|inner:&RetainedError<ck::PdbModelTransitionError>|match inner.as_ref(){ck::PdbModelTransitionError::Structure(value)=>Some(BioStructureError::from_rust(value.clone())),_=>None},
source=Option<crate::error_framework::ErrorObject<BioStructureError>>,methods={
#[getter]fn chain(&self)->Option<&str>{match self.inner.as_ref(){ck::PdbModelTransitionError::NonpolymerSubchainCounterOutsideSourceIntRange{chain}=>Some(chain.as_str()),_=>None}}
#[getter]fn id(&self)->Option<u32>{match self.inner.as_ref(){ck::PdbModelTransitionError::ModelIdOutOfBounds{id,..}=>Some(id.value()),_=>None}}
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbChargeError {
    inner: ck::PdbChargeError,
}
#[pymethods]
impl PdbChargeError {
    #[getter]
    fn digit(&self) -> u8 {
        self.inner.digit()
    }
    #[getter]
    fn sign(&self) -> u8 {
        self.inner.sign()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

value_kind!(PdbAtomFieldErrorKind,PdbAtomFieldError,{TooShort{len},SerialOutsideSourceDefinedRange,ChainNameNotRepresentable{bytes},AtomNameNotRepresentable{bytes},ResidueKey(_),Charge(_)});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbAtomFieldError {
    owner: RetainedError<ck::PdbReaderCause>,
}
impl PdbAtomFieldError {
    fn inner(&self) -> &ck::PdbAtomFieldError {
        let ck::PdbReaderCause::AtomFields(inner) = self.owner.as_ref() else {
            unreachable!("actual native atom-field payload")
        };
        inner
    }
}
#[pymethods]
impl PdbAtomFieldError {
    #[getter]
    fn variant(&self) -> PdbAtomFieldErrorKind {
        self.inner().into()
    }
    #[getter]
    fn len(&self) -> Option<usize> {
        match self.inner() {
            ck::PdbAtomFieldError::TooShort { len } => Some(*len),
            _ => None,
        }
    }
    #[getter]
    fn bytes(&self) -> Option<Source<[u8; 2], [u8; 4]>> {
        match self.inner() {
            ck::PdbAtomFieldError::ChainNameNotRepresentable { bytes } => {
                Some(Source::First(*bytes))
            }
            ck::PdbAtomFieldError::AtomNameNotRepresentable { bytes } => {
                Some(Source::Second(*bytes))
            }
            _ => None,
        }
    }
    #[getter]
    fn residue_key(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match self.inner() {
            ck::PdbAtomFieldError::ResidueKey(value) => PdbResidueKeyError::from_rust(*value)
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
    #[getter]
    fn charge(&self) -> Option<PdbChargeError> {
        match self.inner() {
            ck::PdbAtomFieldError::Charge(inner) => Some(PdbChargeError { inner: *inner }),
            _ => None,
        }
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner())
    }
}
value_kind!(PdbStructureFinalizationErrorKind,PdbStructureFinalizationError,{
Model(_),Connection(_),Remarks(_),Ccd(_),NoModelAfterDefault,SourceIdNotRepresentable{field},Utf8{field,valid_up_to},Structure(_),StagingRowOrder{table,expected,actual}
});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "cosmolkit", frozen)]
pub struct PdbStructureFinalizationError {
    owner: RetainedError<ck::PdbReaderCause>,
}
impl PdbStructureFinalizationError {
    fn inner(&self) -> &ck::PdbStructureFinalizationError {
        let ck::PdbReaderCause::Finalization(inner) = self.owner.as_ref() else {
            unreachable!("actual native finalization payload")
        };
        inner
    }
}
#[pymethods]
impl PdbStructureFinalizationError {
    #[getter]
    fn variant(&self) -> PdbStructureFinalizationErrorKind {
        self.inner().into()
    }
    #[getter]
    fn model(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbModelTransitionError>>> {
        match self.inner() {
            ck::PdbStructureFinalizationError::Model(_) => {
                PdbModelTransitionError::from_ref(self.owner.node().project(|root| {
                    let root = root
                        .downcast_ref::<ck::PdbReaderCause>()
                        .expect("native PDB root");
                    let ck::PdbReaderCause::Finalization(ck::PdbStructureFinalizationError::Model(
                        error,
                    )) = root
                    else {
                        unreachable!("native finalization model payload")
                    };
                    error
                }))
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
    #[getter]
    fn connection(&self) -> Option<PdbConnectionFinalizationError> {
        match self.inner() {
            ck::PdbStructureFinalizationError::Connection(inner) => {
                Some(PdbConnectionFinalizationError {
                    inner: inner.clone(),
                })
            }
            _ => None,
        }
    }
    #[getter]
    fn remarks(&self) -> Option<PdbRemarkFinalizationError> {
        match self.inner() {
            ck::PdbStructureFinalizationError::Remarks(inner) => {
                Some(PdbRemarkFinalizationError { inner: *inner })
            }
            _ => None,
        }
    }
    #[getter]
    fn ccd(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<PdbResidueKeyError>>> {
        match self.inner() {
            ck::PdbStructureFinalizationError::Ccd(inner) => PdbResidueKeyError::from_rust(*inner)
                .value(py)
                .extract()
                .map(Some)
                .map_err(Into::into),
            _ => Ok(None),
        }
    }
    #[getter]
    fn structure(
        &self,
        py: Python<'_>,
    ) -> PyResult<Option<crate::error_framework::ErrorObject<BioStructureError>>> {
        match self.inner() {
            ck::PdbStructureFinalizationError::Structure(inner) => {
                BioStructureError::from_rust(inner.clone())
                    .value(py)
                    .extract()
                    .map(Some)
                    .map_err(Into::into)
            }
            _ => Ok(None),
        }
    }
    #[getter]
    fn field(&self) -> Option<&'static str> {
        match self.inner() {
            ck::PdbStructureFinalizationError::SourceIdNotRepresentable { field }
            | ck::PdbStructureFinalizationError::Utf8 { field, .. } => Some(field),
            _ => None,
        }
    }
    #[getter]
    fn valid_up_to(&self) -> Option<usize> {
        match self.inner() {
            ck::PdbStructureFinalizationError::Utf8 { valid_up_to, .. } => Some(*valid_up_to),
            _ => None,
        }
    }
    #[getter]
    fn table(&self) -> Option<&'static str> {
        match self.inner() {
            ck::PdbStructureFinalizationError::StagingRowOrder { table, .. } => Some(table),
            _ => None,
        }
    }
    #[getter]
    fn expected(&self) -> Option<usize> {
        match self.inner() {
            ck::PdbStructureFinalizationError::StagingRowOrder { expected, .. } => Some(*expected),
            _ => None,
        }
    }
    #[getter]
    fn actual(&self) -> Option<usize> {
        match self.inner() {
            ck::PdbStructureFinalizationError::StagingRowOrder { actual, .. } => Some(*actual),
            _ => None,
        }
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner())
    }
}

error_projection!(retained PdbReaderCause,PdbReaderCauseKind,ck::PdbReaderCause,"bio",{
LineNumberOutsideSourceIntRange,InputFormat(_),AtomFields(_),ModelTransition(_),Remark(_),Conect(_),ModRes(_),Helix(_),Sheet(_),Header(_),DbRef(_),Mtrix(_),Crystal(_),Finalization(_)
},cause=|_:&RetainedError<ck::PdbReaderCause>|None::<PyErr>,methods={
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn input_format(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PdbInputFormatError>>>{match self.inner.as_ref(){ck::PdbReaderCause::InputFormat(value)=>PdbInputFormatError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn atom_fields(&self)->Option<PdbAtomFieldError>{match self.inner.as_ref(){ck::PdbReaderCause::AtomFields(_)=>Some(PdbAtomFieldError{owner:self.inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn model_transition(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PdbModelTransitionError>>>{match self.inner.as_ref(){ck::PdbReaderCause::ModelTransition(_)=>PdbModelTransitionError::from_ref(self.inner.node().project(|root|{let root=root.downcast_ref::<ck::PdbReaderCause>().expect("native PDB cause");let ck::PdbReaderCause::ModelTransition(error)=root else{unreachable!("native model_transition error payload")};error})).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn remark(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PdbRemarkError>>>{match self.inner.as_ref(){ck::PdbReaderCause::Remark(value)=>PdbRemarkError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn conect(&self)->Option<PdbConectError>{match self.inner.as_ref(){ck::PdbReaderCause::Conect(inner)=>Some(PdbConectError{inner:inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn mod_res(&self)->Option<PdbModResError>{match self.inner.as_ref(){ck::PdbReaderCause::ModRes(inner)=>Some(PdbModResError{inner:inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn helix(&self)->Option<PdbHelixError>{match self.inner.as_ref(){ck::PdbReaderCause::Helix(inner)=>Some(PdbHelixError{inner:inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn sheet(&self)->Option<PdbSheetError>{match self.inner.as_ref(){ck::PdbReaderCause::Sheet(inner)=>Some(PdbSheetError{inner:inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn header(&self,py:Python<'_>)->PyResult<Option<crate::error_framework::ErrorObject<PdbHeaderError>>>{match self.inner.as_ref(){ck::PdbReaderCause::Header(value)=>PdbHeaderError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),_=>Ok(None)}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn db_ref(&self)->Option<PdbDbRefError>{match self.inner.as_ref(){ck::PdbReaderCause::DbRef(inner)=>Some(PdbDbRefError{inner:inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn mtrix(&self)->Option<PdbMtrixError>{match self.inner.as_ref(){ck::PdbReaderCause::Mtrix(inner)=>Some(PdbMtrixError{inner:inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn crystal(&self)->Option<PdbCrystalError>{match self.inner.as_ref(){ck::PdbReaderCause::Crystal(inner)=>Some(PdbCrystalError{inner:inner.clone()}),_=>None}}
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn finalization(&self)->Option<PdbStructureFinalizationError>{match self.inner.as_ref(){ck::PdbReaderCause::Finalization(_)=>Some(PdbStructureFinalizationError{owner:self.inner.clone()}),_=>None}}
});

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct BioPdbReadError {
    inner: RetainedError<ck::BioPdbReadError>,
}
}
impl BioPdbReadError {
    pub(crate) fn from_rust(inner: ck::BioPdbReadError) -> PyErr {
        Self::from_retained(RetainedError::new(inner))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Self::from_retained(RetainedError::from_ref(node))
    }
    fn from_retained(inner: RetainedError<ck::BioPdbReadError>) -> PyErr {
        let cause =
            PdbReaderCause::from_ref(inner.node().source().expect("native PDB reader source"));
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => {
                let error = crate::error_framework::from_payload(value);
                error.set_cause(py, Some(cause));
                error
            }
            Err(error) => error,
        })
    }
}
crate::error_framework::error_projection! {methods impl BioPdbReadError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn line_number(&self) -> Option<i32> {
        self.inner.line_number()
    }
    #[getter]
    fn stage(&self) -> BioPdbReadStage {
        self.inner.stage().into()
    }
    #[getter]
    fn record_tag(&self) -> Option<[u8; 4]> {
        self.inner.record_tag()
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<crate::error_framework::ErrorObject<PdbReaderCause>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

pub(crate) fn keywords(
    max_line_length: Option<i32>,
    check_non_ascii: Option<bool>,
    ignore_ter: Option<bool>,
    split_chain_on_ter: Option<bool>,
    skip_remarks: Option<bool>,
) -> crate::params::BioPdbReadParams {
    let defaults = ck::BioPdbReadParams::default();
    crate::params::BioPdbReadParams {
        inner: crate::parameters::NativeField::Owned(ck::BioPdbReadParams {
            max_line_length: max_line_length.unwrap_or(defaults.max_line_length),
            check_non_ascii: check_non_ascii.unwrap_or(defaults.check_non_ascii),
            ignore_ter: ignore_ter.unwrap_or(defaults.ignore_ter),
            split_chain_on_ter: split_chain_on_ter.unwrap_or(defaults.split_chain_on_ter),
            skip_remarks: skip_remarks.unwrap_or(defaults.skip_remarks),
        }),
    }
}
#[pyfunction]
#[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,max_line_length=crate::arguments::Argument::Missing,check_non_ascii=crate::arguments::Argument::Missing,ignore_ter=crate::arguments::Argument::Missing,split_chain_on_ter=crate::arguments::Argument::Missing,skip_remarks=crate::arguments::Argument::Missing))]
pub fn bio_from_pdb(
    py: Python<'_>,
    text: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioPdbReadParams>>>,
    max_line_length: crate::arguments::Argument<Option<i32>>,
    check_non_ascii: crate::arguments::Argument<Option<bool>>,
    ignore_ter: crate::arguments::Argument<Option<bool>>,
    split_chain_on_ter: crate::arguments::Argument<Option<bool>>,
    skip_remarks: crate::arguments::Argument<Option<bool>>,
) -> PyResult<BioStructure> {
    let keyword_form = max_line_length.is_supplied()
        || check_non_ascii.is_supplied()
        || ignore_ter.is_supplied()
        || split_chain_on_ter.is_supplied()
        || skip_remarks.is_supplied();
    crate::arguments::configuration_arguments!(params; max_line_length, check_non_ascii, ignore_ter, split_chain_on_ter, skip_remarks);
    let params = params.as_deref();
    if keyword_form {
        _bio_from_pdb_keywords(
            py,
            text,
            max_line_length,
            check_non_ascii,
            ignore_ter,
            split_chain_on_ter,
            skip_remarks,
        )
    } else {
        _bio_from_pdb_params(py, text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"bio_from_pdb"}"#]
#[pyo3(signature=(text,params=None))]
pub fn _bio_from_pdb_params(
    py: Python<'_>,
    text: &str,
    params: Option<&crate::params::BioPdbReadParams>,
) -> PyResult<BioStructure> {
    let params = pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

    py.detach(|| match params {
        Some(params) => ck::BioStructure::from_pdb_with_params(text, &params),
        None => ck::bio_from_pdb(text),
    })
    .map(|inner| BioStructure { inner })
    .map_err(BioPdbReadError::from_rust)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"bio_from_pdb"}"#]
#[pyo3(signature=(text,*,max_line_length=None,check_non_ascii=None,ignore_ter=None,split_chain_on_ter=None,skip_remarks=None))]
pub fn _bio_from_pdb_keywords(
    py: Python<'_>,
    text: &str,
    max_line_length: Option<i32>,
    check_non_ascii: Option<bool>,
    ignore_ter: Option<bool>,
    split_chain_on_ter: Option<bool>,
    skip_remarks: Option<bool>,
) -> PyResult<BioStructure> {
    let params = keywords(
        max_line_length,
        check_non_ascii,
        ignore_ter,
        split_chain_on_ter,
        skip_remarks,
    );
    _bio_from_pdb_params(py, text, Some(&params))
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::arguments::Argument::{Missing, Value};
    use pyo3::types::PyDict;
    #[test]
    fn pdb_reader_options_conveniences_and_complete_typed_failures_preserve_native_values() {
        Python::initialize();
        Python::attach(|py|->PyResult<()>{
let text="ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nTER\nHETATM    2  O   HOH A   2       1.000   0.000   0.000  1.00 20.00           O  \nHETATM    3 NA    NA B   1       2.000   0.000   0.000  1.00 20.00          NA  \nEND\n";
let default=ck::BioStructure::from_pdb(text).unwrap();
assert_eq!(default.atoms().len(),3);
let expected=default.to_mmcif().unwrap();
let output=py.get_type::<BioStructure>().call_method1("from_pdb",(text,))?.extract::<Py<BioStructure>>()?;
assert_eq!(output.borrow(py).inner.to_mmcif().unwrap(),expected);
assert_eq!(bio_from_pdb(py,text,Missing,Missing,Missing,Missing,Missing,Missing)?.inner.to_mmcif().unwrap(),expected);
for flags in 0..16{
let native=ck::BioPdbReadParams{max_line_length:80,check_non_ascii:flags&1!=0,ignore_ter:flags&2!=0,split_chain_on_ter:flags&4!=0,skip_remarks:flags&8!=0};
let expected=ck::BioStructure::from_pdb_with_params(text,&native).unwrap().to_mmcif().unwrap();
let params=Py::new(py,crate::params::BioPdbReadParams{inner:crate::parameters::NativeField::Owned(native)})?;
let kw=PyDict::new(py);
kw.set_item("max_line_length",native.max_line_length)?;
kw.set_item("check_non_ascii",native.check_non_ascii)?;
kw.set_item("ignore_ter",native.ignore_ter)?;
kw.set_item("split_chain_on_ter",native.split_chain_on_ter)?;
kw.set_item("skip_remarks",native.skip_remarks)?;
for form in 0..3{
let output=match form{0=>py.get_type::<BioStructure>().call_method1("from_pdb",(text,params.clone_ref(py)))?,1=>py.get_type::<BioStructure>().call_method("from_pdb",(text,),Some(&kw))?,_=>py.get_type::<BioStructure>().call_method1("from_pdb_with_params",(text,params.clone_ref(py)))?}.extract::<Py<BioStructure>>()?;
assert_eq!(output.borrow(py).inner.to_mmcif().unwrap(),expected);
}
let native_params=Py::new(py,crate::params::BioPdbReadParams{inner:crate::parameters::NativeField::Owned(native)})?;
let free=bio_from_pdb(py,text,Value(Some(native_params.borrow(py))),Missing,Missing,Missing,Missing,Missing)?;
assert_eq!(free.inner.to_mmcif().unwrap(),expected);
let free=bio_from_pdb(py,text,Missing,Value(Some(80)),Value(Some(native.check_non_ascii)),Value(Some(native.ignore_ter)),Value(Some(native.split_chain_on_ter)),Value(Some(native.skip_remarks)))?;
assert_eq!(free.inner.to_mmcif().unwrap(),expected);
let mixed=py.get_type::<BioStructure>().call_method("from_pdb",(text,params),Some(&kw)).unwrap_err();
assert!(mixed.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
assert!(bio_from_pdb(py,text,Value(Some(native_params.borrow(py))),Value(Some(80)),Missing,Missing,Missing,Missing).err().expect("mixed PDB config rejected").is_instance_of::<pyo3::exceptions::PyTypeError>(py));
}
let kw=PyDict::new(py);kw.set_item("skip_remarks",true)?;
let output=py.get_type::<BioStructure>().call_method("from_pdb",(text,),Some(&kw))?.extract::<Py<BioStructure>>()?;
let native=ck::BioStructure::from_pdb_with_params(text,&ck::BioPdbReadParams{skip_remarks:true,..Default::default()}).unwrap();
assert_eq!(output.borrow(py).inner.to_mmcif().unwrap(),native.to_mmcif().unwrap());
let expected=ck::BioStructure::from_pdb("ATOM").unwrap_err();

let error=py.get_type::<BioStructure>().call_method1("from_pdb",("ATOM",)).unwrap_err();
assert!(error.is_instance(py, <BioPdbReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));assert_eq!(error.value(py).str()?.to_str()?,expected.to_string());
assert_eq!(error.value(py).call_method0("line_number")?.extract::<Option<i32>>()?,expected.line_number());
assert_eq!(error.value(py).call_method0("record_tag")?.extract::<Option<[u8;4]>>()?,expected.record_tag());
assert!(error.value(py).call_method0("stage")?.extract::<BioPdbReadStage>()?==expected.stage().into());
let cause=error.cause(py).expect("real native PDB cause");
assert!(cause.is_instance(py, <PdbReaderCause as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));assert!(cause.cause(py).is_none());
assert!(error.value(py).getattr("source")?.is(cause.value(py)));
let fields=cause.value(py).getattr("atom_fields")?.unbind();
drop(cause);drop(error);
assert!(fields.bind(py).getattr("variant")?.extract::<PdbAtomFieldErrorKind>()?==PdbAtomFieldErrorKind::TooShort);
assert_eq!(fields.bind(py).getattr("len")?.extract::<usize>()?,4);
assert!(fields.bind(py).setattr("len",0).is_err());

let mut bytes=text.lines().next().unwrap().as_bytes().to_vec();
assert!(bytes.len()>=80);bytes[78]=b'x';bytes[79]=b'+';
// Native read_charge deliberately accepts a non-digit first byte as zero.
let accepted=String::from_utf8(bytes.clone()).unwrap();
let expected=ck::BioStructure::from_pdb(&accepted).unwrap().to_mmcif().unwrap();
assert_eq!(bio_from_pdb(py,&accepted,Missing,Missing,Missing,Missing,Missing,Missing)?.inner.to_mmcif().unwrap(),expected);
// Only a digit followed by an invalid sign triggers the actual charge error.
bytes[78]=b'5';bytes[79]=b'x';
let invalid=String::from_utf8(bytes).unwrap();
let native_error=ck::BioStructure::from_pdb(&invalid).unwrap_err();
let native_cause=native_error.source().unwrap().downcast_ref::<ck::PdbReaderCause>().unwrap();
let ck::PdbReaderCause::AtomFields(ck::PdbAtomFieldError::Charge(native_charge))=native_cause else{panic!("actual source charge failure expected")};
let digit=native_charge.digit();let sign=native_charge.sign();
let error=bio_from_pdb(py,&invalid,Missing,Missing,Missing,Missing,Missing,Missing).err().expect("actual source charge failure");
assert_eq!(error.value(py).str()?.to_str()?,native_error.to_string());
let charge=error.cause(py).unwrap().value(py).getattr("atom_fields")?.getattr("charge")?.unbind();
drop(error);drop(invalid);
assert_eq!(charge.bind(py).getattr("digit")?.extract::<u8>()?,digit);
assert_eq!(charge.bind(py).getattr("sign")?.extract::<u8>()?,sign);
assert!(charge.bind(py).setattr("sign",0).is_err());

let source=ck::PdbReaderCause::Finalization(ck::PdbStructureFinalizationError::Model(ck::PdbModelTransitionError::Structure(ck::BioStructureError::RowIndexTooLarge{value:usize::MAX})));
assert!(source.source().is_none());
let parent=PdbReaderCause::from_rust(source);
assert!(parent.cause(py).is_none());
let finalization=parent.value(py).getattr("finalization")?.unbind();
let model=finalization.bind(py).getattr("model")?.unbind();
let native_cause=model.bind(py).getattr("__cause__")?.unbind();
assert!(native_cause.bind(py).is_instance(<BioStructureError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
assert!(model.bind(py).getattr("source")?.is(native_cause.bind(py)));
drop(parent);drop(finalization);drop(model);
assert_eq!(native_cause.bind(py).getattr("value")?.extract::<usize>()?,usize::MAX);

let model=PdbModelTransitionError::from_rust(ck::PdbModelTransitionError::ModelIdOutOfBounds{id:ck::BioModelId::new(u32::MAX),model_count:usize::MAX});
assert!(model.cause(py).is_none());
assert_eq!(model.value(py).getattr("id")?.extract::<u32>()?,u32::MAX);
assert_eq!(model.value(py).getattr("model_count")?.extract::<usize>()?,usize::MAX);
assert!(model.value(py).setattr("id",0).is_err());

Ok(())}).unwrap();
    }
}
