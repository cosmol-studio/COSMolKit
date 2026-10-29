//! Native in-memory drawing and precise SVG/XML failures.
use crate::{
    Molecule,
    errors::child,
    errors::{enum_error, error_projection, exception},
    input,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub struct XmlTextPos(ck::XmlTextPos);
#[wasm_bindgen]
impl XmlTextPos {
    #[wasm_bindgen(getter)]
    pub fn row(&self) -> u32 {
        self.0.row
    }
    #[wasm_bindgen(getter)]
    pub fn col(&self) -> u32 {
        self.0.col
    }
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "number | string")]
    pub type XmlToken;
    #[wasm_bindgen(typescript_type = "XmlParseError")]
    pub type XmlPayload;
}
enum_error!(XmlParseError,XmlParseErrorKind,ck::XmlParseError,{InvalidXmlPrefixUri=>ck::XmlParseError::InvalidXmlPrefixUri(_),UnexpectedXmlUri=>ck::XmlParseError::UnexpectedXmlUri(_),UnexpectedXmlnsUri=>ck::XmlParseError::UnexpectedXmlnsUri(_),InvalidElementNamePrefix=>ck::XmlParseError::InvalidElementNamePrefix(_),DuplicatedNamespace=>ck::XmlParseError::DuplicatedNamespace(_,_),UnknownNamespace=>ck::XmlParseError::UnknownNamespace(_,_),UnexpectedCloseTag=>ck::XmlParseError::UnexpectedCloseTag(_,_,_),UnexpectedEntityCloseTag=>ck::XmlParseError::UnexpectedEntityCloseTag(_),UnknownEntityReference=>ck::XmlParseError::UnknownEntityReference(_,_),MalformedEntityReference=>ck::XmlParseError::MalformedEntityReference(_),EntityReferenceLoop=>ck::XmlParseError::EntityReferenceLoop(_),InvalidAttributeValue=>ck::XmlParseError::InvalidAttributeValue(_),DuplicatedAttribute=>ck::XmlParseError::DuplicatedAttribute(_,_),NoRootNode=>ck::XmlParseError::NoRootNode,UnclosedRootNode=>ck::XmlParseError::UnclosedRootNode,UnexpectedDeclaration=>ck::XmlParseError::UnexpectedDeclaration(_),DtdDetected=>ck::XmlParseError::DtdDetected,NodesLimitReached=>ck::XmlParseError::NodesLimitReached,AttributesLimitReached=>ck::XmlParseError::AttributesLimitReached,NamespacesLimitReached=>ck::XmlParseError::NamespacesLimitReached,InvalidName=>ck::XmlParseError::InvalidName(_),NonXmlChar=>ck::XmlParseError::NonXmlChar(_,_),InvalidChar=>ck::XmlParseError::InvalidChar(_,_,_),InvalidChar2=>ck::XmlParseError::InvalidChar2(_,_,_),InvalidString=>ck::XmlParseError::InvalidString(_,_),InvalidExternalID=>ck::XmlParseError::InvalidExternalID(_),EntityResolver=>ck::XmlParseError::EntityResolver(_,_),InvalidComment=>ck::XmlParseError::InvalidComment(_),InvalidCharacterData=>ck::XmlParseError::InvalidCharacterData(_),UnknownToken=>ck::XmlParseError::UnknownToken(_),UnexpectedEndOfStream=>ck::XmlParseError::UnexpectedEndOfStream},{
#[wasm_bindgen(getter)]pub fn position(&self)->XmlTextPos{XmlTextPos(self.inner.pos())}
#[wasm_bindgen(getter)]pub fn namespace(&self)->Option<String>{match self.inner.as_ref(){ck::XmlParseError::DuplicatedNamespace(value,_)=>Some(value.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn prefix(&self)->Option<String>{match self.inner.as_ref(){ck::XmlParseError::UnknownNamespace(value,_)=>Some(value.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn entity(&self)->Option<String>{match self.inner.as_ref(){ck::XmlParseError::UnknownEntityReference(value,_)=>Some(value.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn attribute(&self)->Option<String>{match self.inner.as_ref(){ck::XmlParseError::DuplicatedAttribute(value,_)=>Some(value.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn message(&self)->Option<String>{match self.inner.as_ref(){ck::XmlParseError::EntityResolver(_,value)=>Some(value.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn character(&self)->Option<String>{match self.inner.as_ref(){ck::XmlParseError::NonXmlChar(v,_)=>Some(v.to_string()),_=>None}}
#[wasm_bindgen(getter)]pub fn expected(&self)->Option<XmlToken>{match self.inner.as_ref(){ck::XmlParseError::UnexpectedCloseTag(v,_,_)=>Some(JsValue::from(v.to_string()).unchecked_into()),ck::XmlParseError::InvalidChar2(v,_,_)=>Some(JsValue::from(v.to_string()).unchecked_into()),ck::XmlParseError::InvalidString(v,_)=>Some(JsValue::from(v.to_string()).unchecked_into()),ck::XmlParseError::InvalidChar(v,_,_)=>Some(JsValue::from(*v).unchecked_into()),_=>None}}
#[wasm_bindgen(getter)]pub fn actual(&self)->Option<XmlToken>{match self.inner.as_ref(){ck::XmlParseError::UnexpectedCloseTag(_,v,_)=>Some(JsValue::from(v.to_string()).unchecked_into()),ck::XmlParseError::InvalidChar(_,v,_)|ck::XmlParseError::InvalidChar2(_,v,_)=>Some(JsValue::from(*v).unchecked_into()),_=>None}}
});
enum_error!(SvgParseError,SvgParseErrorKind,ck::SvgParseError,{NotAnUtf8Str=>ck::SvgParseError::NotAnUtf8Str,SvgzFeatureNotEnabled=>ck::SvgParseError::SvgzFeatureNotEnabled,MalformedGZip=>ck::SvgParseError::MalformedGZip,ElementsLimitReached=>ck::SvgParseError::ElementsLimitReached,InvalidSize=>ck::SvgParseError::InvalidSize,ParsingFailed=>ck::SvgParseError::ParsingFailed(_)},{#[wasm_bindgen(getter,js_name=parsingFailed)]pub fn parsing_failed(&self)->Option<XmlPayload>{match self.inner.as_ref(){ck::SvgParseError::ParsingFailed(_)=>Some(child!(&self.inner,ck::SvgParseError::ParsingFailed,XmlParseError)),_=>None}}});
crate::enums::enum_value!(DoubleLexicalReadErrorKind,"DoubleLexicalReadErrorKind"{Syntax="syntax",TrailingByte="trailing_byte",Overflow="overflow",FloatingEnvironmentNotModeled="floating_environment_not_modeled"});
#[wasm_bindgen]
pub struct SourceNumericStreamState(ck::SourceNumericStreamState);
#[wasm_bindgen]
impl SourceNumericStreamState {
    #[wasm_bindgen(getter)]
    pub fn eof(&self) -> bool {
        self.0.eof
    }
    #[wasm_bindgen(getter)]
    pub fn fail(&self) -> bool {
        self.0.fail
    }
    #[wasm_bindgen(getter)]
    pub fn bad(&self) -> bool {
        self.0.bad
    }
}
error_projection!(DoubleLexicalReadError, ck::DoubleLexicalReadError, {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> DoubleLexicalReadErrorKind {
        self.inner.kind.into()
    }
    #[wasm_bindgen(getter)]
    pub fn consumed(&self) -> usize {
        self.inner.consumed
    }
    #[wasm_bindgen(getter,js_name=valueBits)]
    pub fn value_bits(&self) -> Option<u64> {
        self.inner.value_bits
    }
    #[wasm_bindgen(getter)]
    pub fn state(&self) -> SourceNumericStreamState {
        SourceNumericStreamState(self.inner.state)
    }
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "AllocationError")]
    #[derive(Clone)]
    pub type AllocationCause;
}
enum_error!(UnsignedStreamArrayError,UnsignedStreamArrayErrorKind,ck::UnsignedStreamArrayError,{MissingFirstValue=>ck::UnsignedStreamArrayError::MissingFirstValue,Allocation=>ck::UnsignedStreamArrayError::Allocation(_)},{},source=AllocationCause,cause=|owner:&ReadonlyValue<ck::UnsignedStreamArrayError>|match owner.as_ref(){ck::UnsignedStreamArrayError::Allocation(_)=>Some(child!(owner,ck::UnsignedStreamArrayError::Allocation,crate::external_errors::AllocationError)),_=>None});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "PropertyValueError | TopologyValidationError | CoordinateValidationError | MappingValidationError | KekulizeError | HydrogenError | WedgeError | ValenceError | Coordinate2DError | FromUtf8Error | SvgParseError | PropertyStringError | UnsignedStreamArrayError | DoubleLexicalReadError | IoError | Utf8Error | ParseFloatError | ExternalError"
    )]
    #[derive(Clone)]
    pub type DrawingCause;
}
fn cause(owner: &ReadonlyValue<ck::DrawingError>) -> Option<DrawingCause> {
    match owner.as_ref() {
        ck::DrawingError::Property(_) => Some(child!(
            owner,
            ck::DrawingError::Property,
            crate::model_errors::PropertyValueError
        )),
        ck::DrawingError::Topology(_) => Some(child!(
            owner,
            ck::DrawingError::Topology,
            crate::model_errors::TopologyValidationError
        )),
        ck::DrawingError::Coordinates(_) => Some(child!(
            owner,
            ck::DrawingError::Coordinates,
            crate::model_errors::CoordinateValidationError
        )),
        ck::DrawingError::Mapping(_) => Some(child!(
            owner,
            ck::DrawingError::Mapping,
            crate::model_errors::MappingValidationError
        )),
        ck::DrawingError::Kekulize(_) => Some(child!(
            owner,
            ck::DrawingError::Kekulize,
            crate::kekulize_errors::KekulizeError
        )),
        ck::DrawingError::Hydrogen(_) => Some(child!(
            owner,
            ck::DrawingError::Hydrogen,
            crate::hydrogen_errors::HydrogenError
        )),
        ck::DrawingError::Wedge(_) => Some(child!(
            owner,
            ck::DrawingError::Wedge,
            crate::wedge_errors::WedgeError
        )),
        ck::DrawingError::Valence(_) => Some(child!(
            owner,
            ck::DrawingError::Valence,
            crate::valence_errors::ValenceError
        )),
        ck::DrawingError::CoordinateGeneration(_) => Some(child!(
            owner,
            ck::DrawingError::CoordinateGeneration,
            crate::depict_errors::Coordinate2DError
        )),
        ck::DrawingError::SvgTextProjection(_) => Some(child!(
            owner,
            ck::DrawingError::SvgTextProjection,
            crate::external_errors::FromUtf8Error
        )),
        ck::DrawingError::SvgParse(_) => {
            Some(child!(owner, ck::DrawingError::SvgParse, SvgParseError))
        }
        ck::DrawingError::PropertyString(_) => Some(child!(
            owner,
            ck::DrawingError::PropertyString,
            crate::property_errors::PropertyStringError
        )),
        ck::DrawingError::VariationArray(_) => Some(child!(
            owner,
            ck::DrawingError::VariationArray,
            UnsignedStreamArrayError
        )),
        ck::DrawingError::DataFieldDouble { .. } => Some(
            DoubleLexicalReadError::from_owner(owner.child(|e| match e {
                ck::DrawingError::DataFieldDouble { source, .. } => source,
                _ => unreachable!("immutable lexical source"),
            }))
            .unchecked_into(),
        ),
        ck::DrawingError::PngEncode(_) => Some(
            crate::external_errors::standard_source(owner.child::<dyn std::error::Error>(
                |e| match e {
                    ck::DrawingError::PngEncode(v) => v.as_ref(),
                    _ => unreachable!("immutable encoder source"),
                },
            ))
            .unchecked_into(),
        ),
        ck::DrawingError::StateRows { .. }
        | ck::DrawingError::HydrogenAppend { .. }
        | ck::DrawingError::InvalidDimensions { .. }
        | ck::DrawingError::PixmapAllocation { .. }
        | ck::DrawingError::VariationIndex { .. } => None,
    }
}
enum_error!(DrawingError,DrawingErrorKind,ck::DrawingError,{Property=>ck::DrawingError::Property(_),Topology=>ck::DrawingError::Topology(_),Coordinates=>ck::DrawingError::Coordinates(_),Mapping=>ck::DrawingError::Mapping(_),Kekulize=>ck::DrawingError::Kekulize(_),Hydrogen=>ck::DrawingError::Hydrogen(_),Wedge=>ck::DrawingError::Wedge(_),Valence=>ck::DrawingError::Valence(_),CoordinateGeneration=>ck::DrawingError::CoordinateGeneration(_),SvgTextProjection=>ck::DrawingError::SvgTextProjection(_),SvgParse=>ck::DrawingError::SvgParse(_),PropertyString=>ck::DrawingError::PropertyString(_),VariationArray=>ck::DrawingError::VariationArray(_),PngEncode=>ck::DrawingError::PngEncode(_),StateRows=>ck::DrawingError::StateRows{..},HydrogenAppend=>ck::DrawingError::HydrogenAppend{..},InvalidDimensions=>ck::DrawingError::InvalidDimensions{..},PixmapAllocation=>ck::DrawingError::PixmapAllocation{..},VariationIndex=>ck::DrawingError::VariationIndex{..},DataFieldDouble=>ck::DrawingError::DataFieldDouble{..}},{
#[wasm_bindgen(getter,js_name=field)]pub fn field(&self)->Option<String>{match self.inner.as_ref(){ck::DrawingError::StateRows{field,..}=>Some(field.to_string()),ck::DrawingError::DataFieldDouble{field,..}=>Some(field.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=actual)]pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::DrawingError::StateRows{actual,..}=>Some(*actual),_=>None}}
#[wasm_bindgen(getter,js_name=expected)]pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::DrawingError::StateRows{expected,..}=>Some(*expected),_=>None}}
#[wasm_bindgen(getter,js_name=row)]pub fn row(&self)->Option<usize>{match self.inner.as_ref(){ck::DrawingError::HydrogenAppend{row,..}=>*row,_=>None}}
#[wasm_bindgen(getter,js_name=reason)]pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::DrawingError::HydrogenAppend{reason,..}=>Some(reason.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=width)]pub fn width(&self)->Option<u32>{match self.inner.as_ref(){ck::DrawingError::InvalidDimensions{width,..}=>Some(*width),ck::DrawingError::PixmapAllocation{width,..}=>Some(*width),_=>None}}
#[wasm_bindgen(getter,js_name=height)]pub fn height(&self)->Option<u32>{match self.inner.as_ref(){ck::DrawingError::InvalidDimensions{height,..}=>Some(*height),ck::DrawingError::PixmapAllocation{height,..}=>Some(*height),_=>None}}
#[wasm_bindgen(getter,js_name=index)]pub fn index(&self)->Option<u32>{match self.inner.as_ref(){ck::DrawingError::VariationIndex{index,..}=>Some(*index),_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){ck::DrawingError::VariationIndex{atom_count,..}=>Some(*atom_count),_=>None}}
#[wasm_bindgen(getter,js_name=group)]pub fn group(&self)->Option<usize>{match self.inner.as_ref(){ck::DrawingError::DataFieldDouble{group,..}=>Some(*group),_=>None}}
},source=DrawingCause,cause=cause);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toSvg)]
    pub fn to_svg(
        &self,
        width: Ts<input::U32Input>,
        height: Ts<input::U32Input>,
    ) -> Result<String, JsValue> {
        self.0
            .to_svg(input::read::<u32>(width)?, input::read::<u32>(height)?)
            .map_err(DrawingError::from_rust)
    }
    #[wasm_bindgen(js_name=toPng)]
    pub fn to_png(
        &self,
        width: Ts<input::U32Input>,
        height: Ts<input::U32Input>,
    ) -> Result<Vec<u8>, JsValue> {
        self.0
            .to_png(input::read::<u32>(width)?, input::read::<u32>(height)?)
            .map_err(DrawingError::from_rust)
    }
}
