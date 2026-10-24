//! Exact SMARTS parsing/CX lowering errors; no independently maintained parser.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{AtomPropertyError, BondValueError, MoleculePropertyError};
use crate::model_errors::StereoGroupError;
use crate::parser_errors::ParserCarrierError;
use crate::property_numeric_errors::PropertyUIntReadError;
use crate::property_values::PropertyValueKind;
use crate::query_model_errors::QueryGraphError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(CxStereoGroupMergeError, CxStereoGroupMergeErrorKind, ck::CxStereoGroupMergeError, "search", {
    StereoGroup(_), TrackerPropertyType { actual: PropertyValueKind => |value: &ck::PropertyValueKind| (*value).into() },
    TrackerCount { groups: usize, hashes: usize },
}, cause = |inner: &ck::CxStereoGroupMergeError| match inner {
    ck::CxStereoGroupMergeError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::CxStereoGroupMergeError::TrackerPropertyType { .. } | ck::CxStereoGroupMergeError::TrackerCount { .. } => None,
}, source = Option<crate::error_framework::ErrorObject<StereoGroupError>>);

type CxSource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    Source<
        crate::error_framework::ErrorObject<CxStereoGroupMergeError>,
        Source<
            crate::error_framework::ErrorObject<PropertyUIntReadError>,
            Source<
                crate::error_framework::ErrorObject<BondValueError>,
                Source<
                    crate::error_framework::ErrorObject<AtomPropertyError>,
                    Source<
                        crate::error_framework::ErrorObject<PropertyStringError>,
                        crate::error_framework::ErrorObject<MoleculePropertyError>,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(CxQueryLoweringError, CxQueryLoweringErrorKind, ck::CxQueryLoweringError, "search", {
    StereoGroup(_), StereoGroupMerge(_), Numeric(_), BondProperty(_), AtomProperty(_), Property(_), MoleculeProperty(_),
    AtomIndex { index }, BondIndex { index }, CoordinateCount { actual: usize, expected: usize },
    BondAtomMismatch { atom, bond }, WedgeAtomMismatch { atom, bond },
    UnsupportedRecord { record: &'static str }, InvalidGraph(_),
}, cause = |inner: &ck::CxQueryLoweringError| match inner {
    ck::CxQueryLoweringError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::CxQueryLoweringError::StereoGroupMerge(value) => Some(CxStereoGroupMergeError::from_rust(value.clone())),
    ck::CxQueryLoweringError::Numeric(value) => Some(PropertyUIntReadError::from_rust(value.clone())),
    ck::CxQueryLoweringError::BondProperty(value) => Some(BondValueError::from_rust(value.clone())),
    ck::CxQueryLoweringError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::CxQueryLoweringError::Property(value) => Some(PropertyStringError::from_rust(*value)),
    ck::CxQueryLoweringError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::CxQueryLoweringError::AtomIndex { .. } | ck::CxQueryLoweringError::BondIndex { .. }
    | ck::CxQueryLoweringError::CoordinateCount { .. } | ck::CxQueryLoweringError::BondAtomMismatch { .. }
    | ck::CxQueryLoweringError::WedgeAtomMismatch { .. } | ck::CxQueryLoweringError::UnsupportedRecord { .. }
    | ck::CxQueryLoweringError::InvalidGraph(_) => None,
}, source = Option<CxSource>, methods = {
    #[getter]
    fn index(&self) -> Option<usize> { match self.inner {
        ck::CxQueryLoweringError::AtomIndex { index } | ck::CxQueryLoweringError::BondIndex { index } => Some(index), _ => None,
    } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::CxQueryLoweringError::BondAtomMismatch { atom, .. } | ck::CxQueryLoweringError::WedgeAtomMismatch { atom, .. } => Some(atom), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::CxQueryLoweringError::BondAtomMismatch { bond, .. } | ck::CxQueryLoweringError::WedgeAtomMismatch { bond, .. } => Some(bond), _ => None,
    } }
    #[getter]
    fn message(&self) -> Option<&str> { match &self.inner { ck::CxQueryLoweringError::InvalidGraph(value) => Some(value), _ => None } }
});

pub(crate) type SmartsSource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    Source<
        crate::error_framework::ErrorObject<CxQueryLoweringError>,
        Source<
            crate::error_framework::ErrorObject<QueryGraphError>,
            Source<
                crate::error_framework::ErrorObject<ParserCarrierError>,
                Source<
                    crate::error_framework::ErrorObject<AtomPropertyError>,
                    Source<
                        crate::error_framework::ErrorObject<BondValueError>,
                        Source<
                            crate::error_framework::ErrorObject<MoleculePropertyError>,
                            crate::error_framework::ErrorObject<TemplateAttachmentOrderError>,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(SmartsParseError, SmartsParseErrorKind, ck::SmartsParseError, "search", {
    StereoGroup(_), MissingRecursiveQueryGraph, CxLowering(_), QueryGraph(_), ParserCarrier(_),
    AtomProperty(_), BondProperty(_), MoleculeProperty(_), UnclosedBracket(_),
    UnexpectedCharacter { position, character: u8, context: String }, UnexpectedEnd(_),
    InvalidAtomPrimitive { position, detail: String }, UnclosedParenthesis(_), UnbalancedRingClosure(_),
    SelfRingClosure { ring, atom: usize }, DuplicateRingBond { ring, begin_atom: usize, end_atom: usize },
    CxSmiles(_), Parse(_), UnsupportedFeature(_), TemplateAttachmentRemap { carrier: usize, source },
}, cause = |inner: &ck::SmartsParseError| match inner {
    ck::SmartsParseError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::SmartsParseError::CxLowering(value) => Some(CxQueryLoweringError::from_rust(value.clone())),
    ck::SmartsParseError::QueryGraph(value) => Some(QueryGraphError::from_rust(value.clone())),
    ck::SmartsParseError::ParserCarrier(value) => Some(ParserCarrierError::from_rust(value.clone())),
    ck::SmartsParseError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::SmartsParseError::BondProperty(value) => Some(BondValueError::from_rust(value.clone())),
    ck::SmartsParseError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::SmartsParseError::TemplateAttachmentRemap { source, .. } => Some(TemplateAttachmentOrderError::from_rust(source.clone())),
    ck::SmartsParseError::MissingRecursiveQueryGraph | ck::SmartsParseError::UnclosedBracket(_)
    | ck::SmartsParseError::UnexpectedCharacter { .. } | ck::SmartsParseError::UnexpectedEnd(_)
    | ck::SmartsParseError::InvalidAtomPrimitive { .. } | ck::SmartsParseError::UnclosedParenthesis(_)
    | ck::SmartsParseError::UnbalancedRingClosure(_) | ck::SmartsParseError::SelfRingClosure { .. }
    | ck::SmartsParseError::DuplicateRingBond { .. } | ck::SmartsParseError::CxSmiles(_)
    | ck::SmartsParseError::Parse(_) | ck::SmartsParseError::UnsupportedFeature(_) => None,
}, source = Option<SmartsSource>, methods = {
    #[getter]
    fn position(&self) -> Option<usize> { match self.inner {
        ck::SmartsParseError::UnclosedBracket(value) | ck::SmartsParseError::UnclosedParenthesis(value) => Some(value),
        ck::SmartsParseError::UnexpectedCharacter { position, .. }
        | ck::SmartsParseError::InvalidAtomPrimitive { position, .. } => Some(position), _ => None,
    } }
    #[getter]
    fn ring(&self) -> Option<u32> { match self.inner {
        ck::SmartsParseError::UnbalancedRingClosure(value) => Some(value),
        ck::SmartsParseError::SelfRingClosure { ring, .. } | ck::SmartsParseError::DuplicateRingBond { ring, .. } => Some(ring), _ => None,
    } }
    #[getter]
    fn message(&self) -> Option<&str> { match &self.inner {
        ck::SmartsParseError::UnexpectedEnd(value) | ck::SmartsParseError::CxSmiles(value)
        | ck::SmartsParseError::Parse(value) => Some(value), _ => None,
    } }
    #[getter]
    fn feature(&self) -> Option<&'static str> { match self.inner { ck::SmartsParseError::UnsupportedFeature(value) => Some(value), _ => None } }
});
