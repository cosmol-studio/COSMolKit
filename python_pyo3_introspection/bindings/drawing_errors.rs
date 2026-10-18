//! Actual renderer errors and complete SVG/XML payloads, with no alternate renderer.
#[cfg(test)]
mod tests {
    use super::*;
    use crate::native::Molecule;
    #[test]
    fn complete_renderers_and_writers_preserve_native_bytes_and_rich_error_owners() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("CCO").unwrap();
            let before = source.to_binary().unwrap();
            let mol = Py::new(
                py,
                Molecule {
                    inner: source.clone(),
                },
            )?;
            assert_eq!(
                mol.bind(py)
                    .call_method1("to_svg", (160, 120))?
                    .extract::<String>()?,
                source.to_svg(160, 120).unwrap()
            );
            assert_eq!(
                mol.bind(py)
                    .call_method1("to_png", (160, 120))?
                    .extract::<Vec<u8>>()?,
                source.to_png(160, 120).unwrap()
            );
            let directory = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("artifacts")
                .join(format!("drawing-complete-{}", std::process::id()));
            std::fs::create_dir_all(&directory).unwrap();
            for (name, extension) in [("write_svg", "svg"), ("write_png", "png")] {
                let path = directory.join(format!("image.{extension}"));
                mol.bind(py).call_method1(name, (&path, 160, 120))?;
                let expected = if extension == "svg" {
                    source.to_svg(160, 120).unwrap().into_bytes()
                } else {
                    source.to_png(160, 120).unwrap()
                };
                assert_eq!(std::fs::read(path).unwrap(), expected);
            }
            let native = source.to_svg(0, 120).unwrap_err();
            let error = mol.bind(py).call_method1("to_svg", (0, 120)).unwrap_err();
            assert!(error.is_instance_of::<DrawingError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, native.to_string());
            assert_eq!(error.value(py).getattr("width")?.extract::<u32>()?, 0);
            assert!(error.value(py).setattr("width", 7).is_err());
            let e = mol
                .bind(py)
                .call_method1(
                    "write_png",
                    (directory.join("absent").join("x.png"), 160, 120),
                )
                .unwrap_err();
            assert!(e.is_instance_of::<DrawingWriteError>(py));
            let cause = e.value(py).getattr("source")?;
            assert!(cause.is_instance_of::<IoError>());
            assert!(cause.is(&e.value(py).getattr("__cause__")?));
            drop(e);
            assert!(
                cause
                    .getattr("raw_os_error")?
                    .extract::<Option<i32>>()?
                    .is_some()
            );
            let xml = ck::XmlParseError::UnexpectedCloseTag(
                "expected".into(),
                "actual".into(),
                ck::XmlTextPos { row: 2, col: 3 },
            );
            let e = DrawingError::from_rust(ck::DrawingError::SvgParse(
                ck::SvgParseError::ParsingFailed(xml),
            ));
            let svg = e.value(py).getattr("source")?;
            let payload = svg.getattr("parsing_failed")?;
            drop(e);
            drop(svg);
            assert!(payload.is_instance_of::<XmlParseError>());
            assert!(payload.getattr("__cause__")?.is_none());
            assert_eq!(
                payload
                    .getattr("position")?
                    .getattr("row")?
                    .extract::<u32>()?,
                2
            );
            assert_eq!(
                payload.getattr("expected")?.extract::<String>()?,
                "expected"
            );
            assert_eq!(payload.getattr("actual")?.extract::<String>()?, "actual");
            assert!(payload.setattr("expected", "no").is_err());
            assert!(payload.getattr("position")?.setattr("row", 9).is_err());
            assert_eq!(mol.borrow(py).inner.to_binary().unwrap(), before);
            Ok(())
        })
        .unwrap();
    }
}

use crate::allocation_errors::UnsignedStreamArrayError;
use crate::canonical_errors::PropertyStringError;
use crate::chemical_values::enum_value;
use crate::coordinate_2d_errors::Coordinate2DError;
use crate::error_framework::{Source, error_projection};
use crate::hydrogen_errors::HydrogenError;
use crate::io_source_errors::{ErrorRef, IoError, StandardSource, standard_source};
use crate::kekulize_errors::KekulizeError;
use crate::model_edit_errors::MappingValidationError;
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::property_values::PropertyValueError;
use crate::retained_error::RetainedError;
use crate::smarts_writer_support_errors::WedgeError;
use crate::utf8_errors::FromUtf8Error;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
use std::error::Error;

#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct XmlTextPos {
    inner: ck::XmlTextPos,
}
#[pymethods]
impl XmlTextPos {
    #[getter]
    fn row(&self) -> u32 {
        self.inner.row
    }
    #[getter]
    fn col(&self) -> u32 {
        self.inner.col
    }
}
#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum XmlParseErrorKind {
    InvalidXmlPrefixUri,
    UnexpectedXmlUri,
    UnexpectedXmlnsUri,
    InvalidElementNamePrefix,
    DuplicatedNamespace,
    UnknownNamespace,
    UnexpectedCloseTag,
    UnexpectedEntityCloseTag,
    UnknownEntityReference,
    MalformedEntityReference,
    EntityReferenceLoop,
    InvalidAttributeValue,
    DuplicatedAttribute,
    NoRootNode,
    UnclosedRootNode,
    UnexpectedDeclaration,
    DtdDetected,
    NodesLimitReached,
    AttributesLimitReached,
    NamespacesLimitReached,
    InvalidName,
    NonXmlChar,
    InvalidChar,
    InvalidChar2,
    InvalidString,
    InvalidExternalID,
    EntityResolver,
    InvalidComment,
    InvalidCharacterData,
    UnknownToken,
    UnexpectedEndOfStream,
}
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct XmlParseError {
    inner: RetainedError<ck::XmlParseError>,
}
impl XmlParseError {
    pub(crate) fn from_rust(inner: ck::XmlParseError) -> PyErr {
        Self::from_retained(RetainedError::new(inner))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Self::from_retained(RetainedError::from_ref(node))
    }
    fn from_retained(inner: RetainedError<ck::XmlParseError>) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(v) => PyErr::from_value(v.into_bound(py).into_any()),
            Err(e) => e,
        })
    }
}
#[pymethods]
impl XmlParseError {
    #[getter]
    fn domain(&self) -> &'static str {
        "depict"
    }
    #[getter]
    fn variant(&self) -> XmlParseErrorKind {
        match self.inner.as_ref() {
            ck::XmlParseError::InvalidXmlPrefixUri(_) => XmlParseErrorKind::InvalidXmlPrefixUri,
            ck::XmlParseError::UnexpectedXmlUri(_) => XmlParseErrorKind::UnexpectedXmlUri,
            ck::XmlParseError::UnexpectedXmlnsUri(_) => XmlParseErrorKind::UnexpectedXmlnsUri,
            ck::XmlParseError::InvalidElementNamePrefix(_) => {
                XmlParseErrorKind::InvalidElementNamePrefix
            }
            ck::XmlParseError::DuplicatedNamespace(_, _) => XmlParseErrorKind::DuplicatedNamespace,
            ck::XmlParseError::UnknownNamespace(_, _) => XmlParseErrorKind::UnknownNamespace,
            ck::XmlParseError::UnexpectedCloseTag(_, _, _) => XmlParseErrorKind::UnexpectedCloseTag,
            ck::XmlParseError::UnexpectedEntityCloseTag(_) => {
                XmlParseErrorKind::UnexpectedEntityCloseTag
            }
            ck::XmlParseError::UnknownEntityReference(_, _) => {
                XmlParseErrorKind::UnknownEntityReference
            }
            ck::XmlParseError::MalformedEntityReference(_) => {
                XmlParseErrorKind::MalformedEntityReference
            }
            ck::XmlParseError::EntityReferenceLoop(_) => XmlParseErrorKind::EntityReferenceLoop,
            ck::XmlParseError::InvalidAttributeValue(_) => XmlParseErrorKind::InvalidAttributeValue,
            ck::XmlParseError::DuplicatedAttribute(_, _) => XmlParseErrorKind::DuplicatedAttribute,
            ck::XmlParseError::NoRootNode => XmlParseErrorKind::NoRootNode,
            ck::XmlParseError::UnclosedRootNode => XmlParseErrorKind::UnclosedRootNode,
            ck::XmlParseError::UnexpectedDeclaration(_) => XmlParseErrorKind::UnexpectedDeclaration,
            ck::XmlParseError::DtdDetected => XmlParseErrorKind::DtdDetected,
            ck::XmlParseError::NodesLimitReached => XmlParseErrorKind::NodesLimitReached,
            ck::XmlParseError::AttributesLimitReached => XmlParseErrorKind::AttributesLimitReached,
            ck::XmlParseError::NamespacesLimitReached => XmlParseErrorKind::NamespacesLimitReached,
            ck::XmlParseError::InvalidName(_) => XmlParseErrorKind::InvalidName,
            ck::XmlParseError::NonXmlChar(_, _) => XmlParseErrorKind::NonXmlChar,
            ck::XmlParseError::InvalidChar(_, _, _) => XmlParseErrorKind::InvalidChar,
            ck::XmlParseError::InvalidChar2(_, _, _) => XmlParseErrorKind::InvalidChar2,
            ck::XmlParseError::InvalidString(_, _) => XmlParseErrorKind::InvalidString,
            ck::XmlParseError::InvalidExternalID(_) => XmlParseErrorKind::InvalidExternalID,
            ck::XmlParseError::EntityResolver(_, _) => XmlParseErrorKind::EntityResolver,
            ck::XmlParseError::InvalidComment(_) => XmlParseErrorKind::InvalidComment,
            ck::XmlParseError::InvalidCharacterData(_) => XmlParseErrorKind::InvalidCharacterData,
            ck::XmlParseError::UnknownToken(_) => XmlParseErrorKind::UnknownToken,
            ck::XmlParseError::UnexpectedEndOfStream => XmlParseErrorKind::UnexpectedEndOfStream,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn position(&self) -> XmlTextPos {
        XmlTextPos {
            inner: self.inner.pos(),
        }
    }
    #[getter]
    fn namespace(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::XmlParseError::DuplicatedNamespace(value, _) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn prefix(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::XmlParseError::UnknownNamespace(value, _) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn entity(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::XmlParseError::UnknownEntityReference(value, _) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn attribute(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::XmlParseError::DuplicatedAttribute(value, _) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn message(&self) -> Option<&str> {
        match self.inner.as_ref() {
            ck::XmlParseError::EntityResolver(_, value) => Some(value),
            _ => None,
        }
    }
    #[getter]
    fn character(&self) -> Option<char> {
        match self.inner.as_ref() {
            ck::XmlParseError::NonXmlChar(value, _) => Some(*value),
            _ => None,
        }
    }
    #[getter]
    fn expected(&self) -> Option<Source<u8, String>> {
        match self.inner.as_ref() {
            ck::XmlParseError::UnexpectedCloseTag(value, _, _) => {
                Some(Source::Second(value.to_string()))
            }
            ck::XmlParseError::InvalidChar(value, _, _) => Some(Source::First(*value)),
            ck::XmlParseError::InvalidChar2(value, _, _) => Some(Source::Second(value.to_string())),
            ck::XmlParseError::InvalidString(value, _) => Some(Source::Second(value.to_string())),
            _ => None,
        }
    }
    #[getter]
    fn actual(&self) -> Option<Source<u8, String>> {
        match self.inner.as_ref() {
            ck::XmlParseError::UnexpectedCloseTag(_, value, _) => {
                Some(Source::Second(value.to_string()))
            }
            ck::XmlParseError::InvalidChar(_, value, _) => Some(Source::First(*value)),
            ck::XmlParseError::InvalidChar2(_, value, _) => Some(Source::First(*value)),
            _ => None,
        }
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

error_projection!(retained SvgParseError,SvgParseErrorKind,ck::SvgParseError,"depict",{
NotAnUtf8Str,SvgzFeatureNotEnabled,MalformedGZip,ElementsLimitReached,InvalidSize,ParsingFailed(_),
},cause=|_:&RetainedError<ck::SvgParseError>|None::<PyErr>,methods={
#[getter]#[doc=r#"{"role":"error_payload"}"#]
fn parsing_failed(&self,py:Python<'_>)->PyResult<Option<Py<XmlParseError>>>{match self.inner.as_ref(){
ck::SvgParseError::ParsingFailed(_)=>Ok(Some(XmlParseError::from_ref(self.inner.node().project(|e|match e.downcast_ref::<ck::SvgParseError>().unwrap(){ck::SvgParseError::ParsingFailed(v)=>v,_=>unreachable!()})).value(py).extract()?)),
_=>Ok(None)}}
});
enum_value!(DoubleLexicalReadErrorKind,{Syntax,TrailingByte,Overflow,FloatingEnvironmentNotModeled});
#[doc = r#"{"role":"error_value"}"#]
#[pyclass(module = "fixed_projection._native", frozen)]
pub struct SourceNumericStreamState {
    inner: ck::SourceNumericStreamState,
}
#[pymethods]
impl SourceNumericStreamState {
    #[getter]
    fn eof(&self) -> bool {
        self.inner.eof
    }
    #[getter]
    fn fail(&self) -> bool {
        self.inner.fail
    }
    #[getter]
    fn bad(&self) -> bool {
        self.inner.bad
    }
}
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native",extends=PyValueError)]
pub struct DoubleLexicalReadError {
    inner: RetainedError<ck::DoubleLexicalReadError>,
}
impl DoubleLexicalReadError {
    pub(crate) fn from_rust(inner: ck::DoubleLexicalReadError) -> PyErr {
        Self::from_ref(ErrorRef::new(std::sync::Arc::new(inner)))
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        let inner = RetainedError::from_ref(node);
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(v) => PyErr::from_value(v.into_bound(py).into_any()),
            Err(e) => e,
        })
    }
}
#[pymethods]
impl DoubleLexicalReadError {
    #[getter]
    fn domain(&self) -> &'static str {
        "property"
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> DoubleLexicalReadErrorKind {
        self.inner.kind.into()
    }
    #[getter]
    fn consumed(&self) -> usize {
        self.inner.consumed
    }
    #[getter]
    fn value_bits(&self) -> Option<u64> {
        self.inner.value_bits
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn state(&self) -> SourceNumericStreamState {
        SourceNumericStreamState {
            inner: self.inner.state,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

type DrawingSource = Source<
    Py<PropertyValueError>,
    Source<
        Py<TopologyValidationError>,
        Source<
            Py<CoordinateValidationError>,
            Source<
                Py<MappingValidationError>,
                Source<
                    Py<KekulizeError>,
                    Source<
                        Py<HydrogenError>,
                        Source<
                            Py<WedgeError>,
                            Source<
                                Py<ValenceError>,
                                Source<
                                    Py<Coordinate2DError>,
                                    Source<
                                        Py<FromUtf8Error>,
                                        Source<
                                            Py<SvgParseError>,
                                            Source<
                                                StandardSource,
                                                Source<
                                                    Py<PropertyStringError>,
                                                    Source<
                                                        Py<UnsignedStreamArrayError>,
                                                        Py<DoubleLexicalReadError>,
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
error_projection!(retained DrawingError,DrawingErrorKind,ck::DrawingError,"depict",{
Property(_),Topology(_),Coordinates(_),Mapping(_),Kekulize(_),Hydrogen(_),Wedge(_),Valence(_),CoordinateGeneration(_),
StateRows{field,actual:usize,expected:usize},HydrogenAppend{row:Option<usize>,reason:&'static str},
InvalidDimensions{width,height},SvgTextProjection(_),SvgParse(_),PixmapAllocation{width,height},PngEncode(_),
PropertyString(_),VariationArray(_),VariationIndex{index:u32,atom_count:usize},DataFieldDouble{group:usize,field,source},
},cause=|inner:&RetainedError<ck::DrawingError>|{let node=inner.node().source()?;Some(match inner.as_ref(){
ck::DrawingError::Property(v)=>PropertyValueError::from_rust(v.clone()),
ck::DrawingError::Topology(v)=>TopologyValidationError::from_rust(v.clone()),
ck::DrawingError::Coordinates(v)=>CoordinateValidationError::from_rust(v.clone()),
ck::DrawingError::Mapping(v)=>MappingValidationError::from_rust(v.clone()),
ck::DrawingError::Kekulize(v)=>KekulizeError::from_rust(v.clone()),
ck::DrawingError::Hydrogen(v)=>HydrogenError::from_rust(v.clone()),
ck::DrawingError::Wedge(v)=>WedgeError::from_rust(v.clone()),
ck::DrawingError::Valence(v)=>ValenceError::from_rust(v.clone()),
ck::DrawingError::CoordinateGeneration(_)=>Coordinate2DError::from_ref(node),
ck::DrawingError::SvgTextProjection(v)=>FromUtf8Error::from_rust(v.clone()),
ck::DrawingError::SvgParse(_)=>SvgParseError::from_ref(node),
ck::DrawingError::PngEncode(_)=>standard_source(node),
ck::DrawingError::PropertyString(v)=>PropertyStringError::from_rust(*v),
ck::DrawingError::VariationArray(_)=>UnsignedStreamArrayError::from_ref(node),
ck::DrawingError::DataFieldDouble{source,..}=>DoubleLexicalReadError::from_rust(source.clone()),
ck::DrawingError::StateRows{..}|ck::DrawingError::HydrogenAppend{..}|ck::DrawingError::InvalidDimensions{..}|ck::DrawingError::PixmapAllocation{..}|ck::DrawingError::VariationIndex{..}=>unreachable!("source-free drawing error"),
})},source=Option<DrawingSource>,methods={
#[getter]fn field(&self)->Option<&'static str>{match self.inner.as_ref(){ck::DrawingError::StateRows{field,..}|ck::DrawingError::DataFieldDouble{field,..}=>Some(*field),_=>None}}
#[getter]fn width(&self)->Option<u32>{match self.inner.as_ref(){ck::DrawingError::InvalidDimensions{width,..}|ck::DrawingError::PixmapAllocation{width,..}=>Some(*width),_=>None}}
#[getter]fn height(&self)->Option<u32>{match self.inner.as_ref(){ck::DrawingError::InvalidDimensions{height,..}|ck::DrawingError::PixmapAllocation{height,..}=>Some(*height),_=>None}}
});
error_projection!(retained DrawingWriteError,DrawingWriteErrorKind,ck::DrawingWriteError,"depict",{
Drawing(_),Io{path,source},
},cause=|inner:&RetainedError<ck::DrawingWriteError>|{let node=inner.node().source()?;Some(match inner.as_ref(){ck::DrawingWriteError::Drawing(_)=>DrawingError::from_ref(node),ck::DrawingWriteError::Io{..}=>IoError::from_ref(node)})},source=Option<Source<Py<DrawingError>,Py<IoError>>>,methods={
#[getter]fn path(&self)->Option<std::path::PathBuf>{match self.inner.as_ref(){ck::DrawingWriteError::Io{path,..}=>Some(path.clone()),_=>None}}
});
