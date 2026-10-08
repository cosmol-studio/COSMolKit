//! Canonical search error vocabulary and recursive causes.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct SubstructMatchError {
    inner: ck::SubstructMatchError,
}
#[wasm_bindgen]
impl SubstructMatchError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "search".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        substruct_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(getter,js_name=branch,unchecked_return_type="string | null")]
    pub fn branch(&self) -> JsValue {
        match &self.inner {
            ck::SubstructMatchError::Unsupported { branch, .. } => (*branch).into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=rdkitFunction,unchecked_return_type="string | null")]
    pub fn rdkit_function(&self) -> JsValue {
        match &self.inner {
            ck::SubstructMatchError::Unsupported { rdkit_function, .. } => (*rdkit_function).into(),
            _ => JsValue::NULL,
        }
    }
}
fn substruct_error_kind(source: &ck::SubstructMatchError) -> &'static str {
    use ck::SubstructMatchError as E;
    match source {
        E::FinalCheckMappingLength { .. } => "FinalCheckMappingLength",
        E::FinalCheckMappingIndex { .. } => "FinalCheckMappingIndex",
        E::FinalCheckInvariant { .. } => "FinalCheckInvariant",
        E::FinalCheckBondEndpoint { .. } => "FinalCheckBondEndpoint",
        E::FinalCheckMissingBond { .. } => "FinalCheckMissingBond",
        E::StereoOrder(..) => "StereoOrder",
        E::PropertyInteger { .. } => "PropertyInteger",
        E::Unsupported { .. } => "Unsupported",
        E::PeriodicTable(..) => "PeriodicTable",
        E::PropertyString(..) => "PropertyString",
        E::QueryContext(..) => "QueryContext",
    }
}
pub(crate) fn substruct_error(source: &ck::SubstructMatchError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SubstructMatchError");
    set(&error, "domain", "search".into())?;
    set(&error, "kind", substruct_error_kind(source).into())?;
    set(
        &error,
        "cause",
        source.source().map_or(Ok(JsValue::NULL), source_error)?,
    )?;
    set(
        &error,
        "detail",
        SubstructMatchError {
            inner: source.clone(),
        }
        .into(),
    )?;
    if let ck::SubstructMatchError::Unsupported {
        branch,
        rdkit_function,
    } = source
    {
        set(&error, "branch", (*branch).into())?;
        set(&error, "rdkitFunction", (*rdkit_function).into())?;
    }
    Ok(error.into())
}
#[wasm_bindgen]
pub struct MatchError {
    inner: ck::MatchError,
}
#[wasm_bindgen]
impl MatchError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "search".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        match_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn match_error_kind(source: &ck::MatchError) -> &'static str {
    use ck::MatchError as E;
    match source {
        E::InvalidQuery(..) => "InvalidQuery",
        E::InvalidTarget(..) => "InvalidTarget",
        E::UnsupportedAtomPredicate(..) => "UnsupportedAtomPredicate",
        E::UnsupportedBondPredicate(..) => "UnsupportedBondPredicate",
        E::Substruct(..) => "Substruct",
    }
}
pub(crate) fn match_error(source: &ck::MatchError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("MatchError");
    set(&error, "domain", "search".into())?;
    set(&error, "kind", match_error_kind(source).into())?;
    set(
        &error,
        "cause",
        source.source().map_or(Ok(JsValue::NULL), source_error)?,
    )?;
    set(
        &error,
        "detail",
        MatchError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(error.into())
}
#[wasm_bindgen]
pub struct QueryCompileError {
    inner: ck::QueryCompileError,
}
#[wasm_bindgen]
impl QueryCompileError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "search".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        compile_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn compile_error_kind(source: &ck::QueryCompileError) -> &'static str {
    use ck::QueryCompileError as E;
    match source {
        E::InvalidGraph(..) => "InvalidGraph",
    }
}
pub(crate) fn compile_error(source: &ck::QueryCompileError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("QueryCompileError");
    set(&error, "domain", "search".into())?;
    set(&error, "kind", compile_error_kind(source).into())?;
    set(
        &error,
        "cause",
        source.source().map_or(Ok(JsValue::NULL), source_error)?,
    )?;
    set(
        &error,
        "detail",
        QueryCompileError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(error.into())
}
#[wasm_bindgen]
pub struct SmartsWriteError {
    inner: ck::SmartsWriteError,
}
#[wasm_bindgen]
impl SmartsWriteError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "search".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        write_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn write_error_kind(source: &ck::SmartsWriteError) -> &'static str {
    use ck::SmartsWriteError as E;
    match source {
        E::Traversal(..) => "Traversal",
        E::CxCoordinates(..) => "CxCoordinates",
        E::CxRingInfo(..) => "CxRingInfo",
        E::CxRingAtomOrderIndex { .. } => "CxRingAtomOrderIndex",
        E::CxRingStereoReferenceMissing { .. } => "CxRingStereoReferenceMissing",
        E::CxWedge(..) => "CxWedge",
        E::CxBondConfigAtropMissingCarriers { .. } => "CxBondConfigAtropMissingCarriers",
        E::CxAtomPropertyOutput(..) => "CxAtomPropertyOutput",
        E::CxMissingConformer => "CxMissingConformer",
        E::CxCoordinateSource(..) => "CxCoordinateSource",
        E::CxCoordinateOutput(..) => "CxCoordinateOutput",
        E::CxSourceBondOutOfRange { .. } => "CxSourceBondOutOfRange",
        E::CxMoleculePropertyUInt { .. } => "CxMoleculePropertyUInt",
        E::CxSgroupVectorCast { .. } => "CxSgroupVectorCast",
        E::CxSgroupPropertyWrite { .. } => "CxSgroupPropertyWrite",
        E::CxStereoGroup(..) => "CxStereoGroup",
        E::CxSourceAtomOutOfRange { .. } => "CxSourceAtomOutOfRange",
        E::MoleculePropertyWrite(..) => "MoleculePropertyWrite",
        E::SourceAtomCount { .. } => "SourceAtomCount",
        E::AtomPropertyWrite { .. } => "AtomPropertyWrite",
        E::CanonicalTraversal(..) => "CanonicalTraversal",
        E::UnwritableBondQuery { .. } => "UnwritableBondQuery",
        E::SourceAtomToLeftIndex { .. } => "SourceAtomToLeftIndex",
        E::SourceBondBeginIndex { .. } => "SourceBondBeginIndex",
        E::AtomMapInt { .. } => "AtomMapInt",
        E::AtomTypeAtomicNumber { .. } => "AtomTypeAtomicNumber",
        E::ChargeMagnitudeOverflow { .. } => "ChargeMagnitudeOverflow",
        E::PropertyValue(..) => "PropertyValue",
        E::CxRequiredProperty { .. } => "CxRequiredProperty",
        E::CxPropertyList { .. } => "CxPropertyList",
        E::CxCoordinateSelectionArity { .. } => "CxCoordinateSelectionArity",
        E::CxMissingOutputOrder { .. } => "CxMissingOutputOrder",
        E::CxOutputOrderPropertyType { .. } => "CxOutputOrderPropertyType",
        E::CxAtomPropertyKind { .. } => "CxAtomPropertyKind",
        E::CxAtomPropertyUInt { .. } => "CxAtomPropertyUInt",
        E::CxAtomPropertyWrite { .. } => "CxAtomPropertyWrite",
        E::CxCoordinateStorage { .. } => "CxCoordinateStorage",
        E::CxCoordinateSelection { .. } => "CxCoordinateSelection",
        E::CxComposition(..) => "CxComposition",
        E::CxOutputOrder { .. } => "CxOutputOrder",
        E::CxRowCount { .. } => "CxRowCount",
        E::CxBondPropertyUInt { .. } => "CxBondPropertyUInt",
        E::CxAtomPropertyInt { .. } => "CxAtomPropertyInt",
        E::CxSgroupPropertyUInt { .. } => "CxSgroupPropertyUInt",
        E::InvalidPropertyKind { .. } => "InvalidPropertyKind",
        E::Property(..) => "Property",
        E::Valence(..) => "Valence",
        E::InvalidGraph(..) => "InvalidGraph",
        E::QueryGraphTraversalUnsupported { .. } => "QueryGraphTraversalUnsupported",
        E::OrAboveAndBelowAnd => "OrAboveAndBelowAnd",
        E::UnknownCombination { .. } => "UnknownCombination",
        E::MissingRecursiveQueryMolecule => "MissingRecursiveQueryMolecule",
        E::SourceBondDirection { .. } => "SourceBondDirection",
        E::UnsupportedBondQuery { .. } => "UnsupportedBondQuery",
        E::UnsupportedAtomQuery { .. } => "UnsupportedAtomQuery",
        E::CompositeChildCount { .. } => "CompositeChildCount",
        E::XorComposite => "XorComposite",
        E::QueryGraphCxExtensionsUnsupported { .. } => "QueryGraphCxExtensionsUnsupported",
        E::RootedAtomOutOfRange { .. } => "RootedAtomOutOfRange",
        E::EmptyAtomSelection => "EmptyAtomSelection",
        E::EmptyBondSelection => "EmptyBondSelection",
        E::FragmentAtomOutOfRange { .. } => "FragmentAtomOutOfRange",
        E::FragmentBondOutOfRange { .. } => "FragmentBondOutOfRange",
        E::BondAtomNotEndpoint { .. } => "BondAtomNotEndpoint",
    }
}
pub(crate) fn write_error(source: &ck::SmartsWriteError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SmartsWriteError");
    set(&error, "domain", "search".into())?;
    set(&error, "kind", write_error_kind(source).into())?;
    set(
        &error,
        "cause",
        source.source().map_or(Ok(JsValue::NULL), source_error)?,
    )?;
    set(
        &error,
        "detail",
        SmartsWriteError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(error.into())
}
