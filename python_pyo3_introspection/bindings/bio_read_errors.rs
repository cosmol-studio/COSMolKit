//! Exact generic BIO reader causes; parser state never crosses the facade.
use crate::bio_errors::BioStructureError;
use crate::bio_mmcif_reader::BioMmcifReadError;
use crate::bio_mmcif_writer::CifReadError;
use crate::bio_pdb_reader::BioPdbReadError;
use crate::bio_values::{BioCoordinateFormat, BioCoordinateFormatInput, BioReadParams};
use crate::chemical_values::enum_value;
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::{ErrorRef, IoError};
use crate::protein_values::ProteinProjectionError;
use crate::retained_error::RetainedError;
use crate::utf8_errors::FromUtf8Error;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;
use std::path::PathBuf;

pub(crate) fn keywords(
    format: Option<BioCoordinateFormatInput>,
    source_name: Option<String>,
) -> PyResult<BioReadParams> {
    let defaults = ck::BioReadParams::default();
    Ok(BioReadParams {
        inner: ck::BioReadParams {
            format: match format {
                Some(format) => format.into_rust()?,
                None => defaults.format,
            },
            source_name: source_name.unwrap_or(defaults.source_name),
        },
    })
}

fn project(node: ErrorRef) -> PyErr {
    macro_rules! retained {($($name:ident),+)=>{$(if node.error().is::<ck::$name>() {return $name::from_ref(node);})+};}
    retained!(
        BioReadError,
        BioPdbReadError,
        BioMmcifReadError,
        MmjsonReadError,
        MmjsonStructureError,
        MmjsonValueError,
        ChemCompModelError,
        ChemCompAtomError
    );
    if node.error().is::<std::io::Error>() {
        return IoError::from_ref(node);
    }
    if let Some(value) = node.error().downcast_ref::<ck::JsonError>() {
        return JsonError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::CifReadError>() {
        return CifReadError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::BioStructureError>() {
        return BioStructureError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<ck::ProteinProjectionError>() {
        return ProteinProjectionError::from_rust(value.clone());
    }
    if let Some(value) = node.error().downcast_ref::<std::string::FromUtf8Error>() {
        return FromUtf8Error::from_rust(value.clone());
    }
    unreachable!("unmodeled actual BIO reader cause")
}

enum_value!(JsonErrorKind, {
    NoError, OutOfMemory, UnexpectedEnd, MissingRootElement, BadRoot, ExpectedComma,
    MissingObjectKey, ExpectedColon, ExpectedEndOfInput, UnexpectedComma, ExpectedValue,
    ExpectedNull, ExpectedFalse, ExpectedTrue, InvalidNumber, MissingExponent,
    IllegalCodepoint, InvalidUnicodeEscape, UnexpectedEndOfUtf16, ExpectedU,
    InvalidUtf16TrailSurrogate, UnknownEscape, InvalidUtf8
});
// Python keywords need binding-level identifier escaping, not stub rewrites.
#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum MmjsonUnexpectedType {
    Null,
    False_,
    True_,
    Array,
    Object,
}
impl From<ck::MmjsonUnexpectedType> for MmjsonUnexpectedType {
    fn from(value: ck::MmjsonUnexpectedType) -> Self {
        match value {
            ck::MmjsonUnexpectedType::Null => Self::Null,
            ck::MmjsonUnexpectedType::False => Self::False_,
            ck::MmjsonUnexpectedType::True => Self::True_,
            ck::MmjsonUnexpectedType::Array => Self::Array,
            ck::MmjsonUnexpectedType::Object => Self::Object,
        }
    }
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct JsonError {
    inner: ck::JsonError,
}
impl JsonError {
    fn from_rust(inner: ck::JsonError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}
#[pymethods]
impl JsonError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> JsonErrorKind {
        self.inner.kind().into()
    }
    #[getter]
    fn byte_offset(&self) -> usize {
        self.inner.byte_offset()
    }
    #[getter]
    fn line(&self) -> usize {
        self.inner.line()
    }
    #[getter]
    fn column(&self) -> usize {
        self.inner.column()
    }
    #[getter]
    fn argument(&self) -> Option<i32> {
        self.inner.argument()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

error_projection!(retained MmjsonValueError, MmjsonValueErrorKind, ck::MmjsonValueError, "bio", {
    UnexpectedType(_), InvalidArrayElement(_), MissingArrayFirstElement,
    InvalidSpan { start: usize, end: usize }, InvalidUtf8 { valid_up_to: usize },
}, cause=|_:&RetainedError<ck::MmjsonValueError>| None::<PyErr>, methods={
    #[getter] #[doc = r#"{"role":"error_payload"}"#]
    fn value(&self) -> Option<MmjsonUnexpectedType> {
        match self.inner.as_ref() {
            ck::MmjsonValueError::UnexpectedType(value) | ck::MmjsonValueError::InvalidArrayElement(value) => Some((*value).into()),
            _ => None,
        }
    }
});
error_projection!(retained MmjsonStructureError, MmjsonStructureErrorKind, ck::MmjsonStructureError, "bio", {
    RootNotObject, BlockKeyPrefix{block}, BlockNotObject{block}, InvalidCategory{block,category:usize},
    ExpectedArray{column,kind:&'static str}, ArrayLength{column,expected:usize,actual:usize},
    RowCountOverflow, MissingArenaValue, Text(_),
}, cause=|inner:&RetainedError<ck::MmjsonStructureError>|inner.node().source().map(project),
source=Option<Py<MmjsonValueError>>, methods={
    #[getter] fn block(&self)->Option<usize> {match self.inner.as_ref(){
        ck::MmjsonStructureError::BlockKeyPrefix{block}|ck::MmjsonStructureError::BlockNotObject{block}
        |ck::MmjsonStructureError::InvalidCategory{block,..}=>Some(*block),_=>None}}
    #[getter] fn column(&self)->Option<usize> {match self.inner.as_ref(){
        ck::MmjsonStructureError::ExpectedArray{column,..}|ck::MmjsonStructureError::ArrayLength{column,..}=>Some(*column),_=>None}}
});
error_projection!(retained MmjsonReadError, MmjsonReadErrorKind, ck::MmjsonReadError, "bio", {
    Parse{source_name,error}, Structure{source_name,error},
}, cause=|inner:&RetainedError<ck::MmjsonReadError>|inner.node().source().map(project),
source=Option<Source<Py<JsonError>,Py<MmjsonStructureError>>>, methods={
    #[getter] fn source_name(&self)->&str {match self.inner.as_ref(){
        ck::MmjsonReadError::Parse{source_name,..}|ck::MmjsonReadError::Structure{source_name,..}=>source_name}}
    #[getter] #[doc = r#"{"role":"error_payload"}"#]
    fn error(slf:Bound<'_,Self>)->PyResult<Source<Py<JsonError>,Py<MmjsonStructureError>>>{
        Ok(slf.getattr("__cause__")?.extract()?)
    }
});
error_projection!(retained ChemCompAtomError, ChemCompAtomErrorKind, ck::ChemCompAtomError,"bio", {
    MissingValue{row,column:usize},NameOutsideCanonicalBoundary{row,value},
    ChargeOutsideSourceDefinedRange{row,value},Number(_),
},cause=|inner:&RetainedError<ck::ChemCompAtomError>|inner.node().source().map(project),
source=Option<Py<CifReadError>>,methods={
    #[getter] fn row(&self)->Option<usize>{match self.inner.as_ref(){
        ck::ChemCompAtomError::MissingValue{row,..}|ck::ChemCompAtomError::NameOutsideCanonicalBoundary{row,..}
        |ck::ChemCompAtomError::ChargeOutsideSourceDefinedRange{row,..}=>Some(*row),_=>None}}
    #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){
        ck::ChemCompAtomError::NameOutsideCanonicalBoundary{value,..}|ck::ChemCompAtomError::ChargeOutsideSourceDefinedRange{value,..}=>Some(value),_=>None}}
});
error_projection!(retained ChemCompModelError,ChemCompModelErrorKind,ck::ChemCompModelError,"bio",{
    NotChemComp,Table(_),Atom(_),Hierarchy(_),ResidueName(_),
},cause=|inner:&RetainedError<ck::ChemCompModelError>|inner.node().source().map(project),
source=Option<Source<Py<CifReadError>,Source<Py<ChemCompAtomError>,Py<BioStructureError>>>>,methods={
    #[getter] fn value(&self)->Option<&str>{match self.inner.as_ref(){ck::ChemCompModelError::ResidueName(value)=>Some(value),_=>None}}
});

type BioReadSource = Source<
    Py<IoError>,
    Source<
        Py<FromUtf8Error>,
        Source<
            Py<BioPdbReadError>,
            Source<
                Py<CifReadError>,
                Source<Py<BioMmcifReadError>, Source<Py<MmjsonReadError>, Py<ChemCompModelError>>>,
            >,
        >,
    >,
>;
error_projection!(retained BioReadError,BioReadErrorKind,ck::BioReadError,"bio",{
    Io{path,source},Utf8{path,source},Pdb(_),Cif(_),Mmcif(_),Mmjson(_),ChemComp(_),
    WrongFormat{source_name:String,format},UnknownFileFormat(_),
},cause=|inner:&RetainedError<ck::BioReadError>|inner.node().source().map(project),
source=Option<BioReadSource>,methods={
    #[getter] fn path(&self)->Option<PathBuf>{match self.inner.as_ref(){
        ck::BioReadError::Io{path,..}|ck::BioReadError::Utf8{path,..}|ck::BioReadError::UnknownFileFormat(path)=>Some(path.clone()),_=>None}}
    #[getter] #[doc = r#"{"role":"error_payload"}"#]
    fn format(&self)->Option<BioCoordinateFormat>{match self.inner.as_ref(){ck::BioReadError::WrongFormat{format,..}=>Some((*format).into()),_=>None}}
});
error_projection!(retained ProteinReadError,ProteinReadErrorKind,ck::ProteinReadError,"bio",{
    Structure(_),Pdb(_),Mmcif(_),Projection(_),
},cause=|inner:&RetainedError<ck::ProteinReadError>|inner.node().source().map(project),
source=Option<Source<Py<BioReadError>,Source<Py<BioPdbReadError>,Source<Py<BioMmcifReadError>,Py<ProteinProjectionError>>>>>);
