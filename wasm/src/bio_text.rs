//! Automatic BIO text-format dispatch retains its concrete native causes.
use crate::{
    bio_hierarchy::{BioCoordinateFormat, BioStructure},
    bio_operations::ProteinProjectionError,
    bio_read_errors::{BioMmcifReadError, BioStructureError, CifReadError},
    errors::{enum_error, error_projection, exception},
    external_errors::{FromUtf8Error, IoError},
    input,
    pdb_read_errors::BioPdbReadError,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
crate::enums::enum_value!(JsonErrorKind,"JsonErrorKind"{NoError="no_error",OutOfMemory="out_of_memory",UnexpectedEnd="unexpected_end",MissingRootElement="missing_root_element",BadRoot="bad_root",ExpectedComma="expected_comma",MissingObjectKey="missing_object_key",ExpectedColon="expected_colon",ExpectedEndOfInput="expected_end_of_input",UnexpectedComma="unexpected_comma",ExpectedValue="expected_value",ExpectedNull="expected_null",ExpectedFalse="expected_false",ExpectedTrue="expected_true",InvalidNumber="invalid_number",MissingExponent="missing_exponent",IllegalCodepoint="illegal_codepoint",InvalidUnicodeEscape="invalid_unicode_escape",UnexpectedEndOfUtf16="unexpected_end_of_utf16",ExpectedU="expected_u",InvalidUtf16TrailSurrogate="invalid_utf16_trail_surrogate",UnknownEscape="unknown_escape",InvalidUtf8="invalid_utf8"});
crate::enums::enum_value!(MmjsonUnexpectedType,"MmjsonUnexpectedType"{Null="null",False="false",True="true",Array="array",Object="object"});
error_projection!(JsonError, ck::JsonError, {
    #[wasm_bindgen(getter,js_name=kind)]
    pub fn kind(&self) -> JsonErrorKind {
        self.inner.kind().into()
    }
    #[wasm_bindgen(getter,js_name=byteOffset)]
    pub fn byte_offset(&self) -> usize {
        self.inner.byte_offset()
    }
    #[wasm_bindgen(getter,js_name=line)]
    pub fn line(&self) -> usize {
        self.inner.line()
    }
    #[wasm_bindgen(getter,js_name=column)]
    pub fn column(&self) -> usize {
        self.inner.column()
    }
    #[wasm_bindgen(getter,js_name=argument)]
    pub fn argument(&self) -> Option<i32> {
        self.inner.argument()
    }
});
enum_error!(MmjsonValueError,MmjsonValueErrorKind,ck::MmjsonValueError,{UnexpectedType=>ck::MmjsonValueError::UnexpectedType(_),InvalidArrayElement=>ck::MmjsonValueError::InvalidArrayElement(_),MissingArrayFirstElement=>ck::MmjsonValueError::MissingArrayFirstElement,InvalidSpan=>ck::MmjsonValueError::InvalidSpan{..},InvalidUtf8=>ck::MmjsonValueError::InvalidUtf8{..}},{
#[wasm_bindgen(getter,js_name=start)]pub fn start(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonValueError::InvalidSpan{start,..}=>Some(*start),_=>None}}
#[wasm_bindgen(getter,js_name=end)]pub fn end(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonValueError::InvalidSpan{end,..}=>Some(*end),_=>None}}
#[wasm_bindgen(getter,js_name=validUpTo)]pub fn valid_up_to(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonValueError::InvalidUtf8{valid_up_to,..}=>Some(*valid_up_to),_=>None}}
#[wasm_bindgen(getter)]pub fn value(&self)->Option<MmjsonUnexpectedType>{match self.inner.as_ref(){ck::MmjsonValueError::UnexpectedType(v)|ck::MmjsonValueError::InvalidArrayElement(v)=>Some((*v).into()),_=>None}}
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "MmjsonValueError")]
    #[derive(Clone)]
    pub type MmjsonStructureErrorCause;
}
enum_error!(MmjsonStructureError,MmjsonStructureErrorKind,ck::MmjsonStructureError,{RootNotObject=>ck::MmjsonStructureError::RootNotObject,BlockKeyPrefix=>ck::MmjsonStructureError::BlockKeyPrefix{..},BlockNotObject=>ck::MmjsonStructureError::BlockNotObject{..},InvalidCategory=>ck::MmjsonStructureError::InvalidCategory{..},ExpectedArray=>ck::MmjsonStructureError::ExpectedArray{..},ArrayLength=>ck::MmjsonStructureError::ArrayLength{..},RowCountOverflow=>ck::MmjsonStructureError::RowCountOverflow,MissingArenaValue=>ck::MmjsonStructureError::MissingArenaValue,Text=>ck::MmjsonStructureError::Text(_)},{
#[wasm_bindgen(getter,js_name=block)]pub fn block(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonStructureError::BlockKeyPrefix{block,..}|ck::MmjsonStructureError::BlockNotObject{block,..}|ck::MmjsonStructureError::InvalidCategory{block,..}=>Some(*block),_=>None}}
#[wasm_bindgen(getter,js_name=category)]pub fn category(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonStructureError::InvalidCategory{category,..}=>Some(*category),_=>None}}
#[wasm_bindgen(getter,js_name=column)]pub fn column(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonStructureError::ExpectedArray{column,..}|ck::MmjsonStructureError::ArrayLength{column,..}=>Some(*column),_=>None}}
#[wasm_bindgen(getter,js_name=kind)]pub fn kind(&self)->Option<String>{match self.inner.as_ref(){ck::MmjsonStructureError::ExpectedArray{kind,..}=>Some(kind.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=expected)]pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonStructureError::ArrayLength{expected,..}=>Some(*expected),_=>None}}
#[wasm_bindgen(getter,js_name=actual)]pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::MmjsonStructureError::ArrayLength{actual,..}=>Some(*actual),_=>None}}
},source=MmjsonStructureErrorCause,cause=|owner:&ReadonlyValue<ck::MmjsonStructureError>|Some(match owner.as_ref(){ck::MmjsonStructureError::Text(_)=>MmjsonValueError::from_owner(owner.child(|e|match e{ck::MmjsonStructureError::Text(v)=>v,_=>unreachable!("immutable native text source")})),_=>return None}.unchecked_into()));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError")]
    #[derive(Clone)]
    pub type ChemCompAtomErrorCause;
}
enum_error!(ChemCompAtomError,ChemCompAtomErrorKind,ck::ChemCompAtomError,{MissingValue=>ck::ChemCompAtomError::MissingValue{..},NameOutsideCanonicalBoundary=>ck::ChemCompAtomError::NameOutsideCanonicalBoundary{..},ChargeOutsideSourceDefinedRange=>ck::ChemCompAtomError::ChargeOutsideSourceDefinedRange{..},Number=>ck::ChemCompAtomError::Number(_)},{
#[wasm_bindgen(getter,js_name=row)]pub fn row(&self)->Option<usize>{match self.inner.as_ref(){ck::ChemCompAtomError::MissingValue{row,..}|ck::ChemCompAtomError::NameOutsideCanonicalBoundary{row,..}|ck::ChemCompAtomError::ChargeOutsideSourceDefinedRange{row,..}=>Some(*row),_=>None}}
#[wasm_bindgen(getter,js_name=column)]pub fn column(&self)->Option<usize>{match self.inner.as_ref(){ck::ChemCompAtomError::MissingValue{column,..}=>Some(*column),_=>None}}
#[wasm_bindgen(getter,js_name=value)]pub fn value(&self)->Option<String>{match self.inner.as_ref(){ck::ChemCompAtomError::NameOutsideCanonicalBoundary{value,..}|ck::ChemCompAtomError::ChargeOutsideSourceDefinedRange{value,..}=>Some(value.to_string()),_=>None}}
},source=ChemCompAtomErrorCause,cause=|owner:&ReadonlyValue<ck::ChemCompAtomError>|Some(match owner.as_ref(){ck::ChemCompAtomError::Number(_)=>CifReadError::from_owner(owner.child(|e|match e{ck::ChemCompAtomError::Number(v)=>v,_=>unreachable!("immutable native text source")})),_=>return None}.unchecked_into()));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "CifReadError | ChemCompAtomError | BioStructureError")]
    #[derive(Clone)]
    pub type ChemCompModelErrorCause;
}
enum_error!(ChemCompModelError,ChemCompModelErrorKind,ck::ChemCompModelError,{NotChemComp=>ck::ChemCompModelError::NotChemComp,Table=>ck::ChemCompModelError::Table(_),Atom=>ck::ChemCompModelError::Atom(_),Hierarchy=>ck::ChemCompModelError::Hierarchy(_),ResidueName=>ck::ChemCompModelError::ResidueName(_)},{
#[wasm_bindgen(getter,js_name=residueName)]pub fn residue_name(&self)->Option<String>{match self.inner.as_ref(){ck::ChemCompModelError::ResidueName(v)=>Some(v.clone()),_=>None}}
},source=ChemCompModelErrorCause,cause=|owner:&ReadonlyValue<ck::ChemCompModelError>|Some(match owner.as_ref(){ck::ChemCompModelError::Table(_)=>CifReadError::from_owner(owner.child(|e|match e{ck::ChemCompModelError::Table(v)=>v,_=>unreachable!("immutable native text source")})),ck::ChemCompModelError::Atom(_)=>ChemCompAtomError::from_owner(owner.child(|e|match e{ck::ChemCompModelError::Atom(v)=>v,_=>unreachable!("immutable native text source")})),ck::ChemCompModelError::Hierarchy(_)=>BioStructureError::from_owner(owner.child(|e|match e{ck::ChemCompModelError::Hierarchy(v)=>v,_=>unreachable!("immutable native text source")})),_=>return None}.unchecked_into()));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "BioReadError | BioPdbReadError | BioMmcifReadError | ProteinProjectionError"
    )]
    #[derive(Clone)]
    pub type ProteinReadErrorCause;
}
enum_error!(ProteinReadError,ProteinReadErrorKind,ck::ProteinReadError,{Structure=>ck::ProteinReadError::Structure(_),Pdb=>ck::ProteinReadError::Pdb(_),Mmcif=>ck::ProteinReadError::Mmcif(_),Projection=>ck::ProteinReadError::Projection(_)},{
},source=ProteinReadErrorCause,cause=|owner:&ReadonlyValue<ck::ProteinReadError>|Some(match owner.as_ref(){ck::ProteinReadError::Structure(_)=>BioReadError::from_owner(owner.child(|e|match e{ck::ProteinReadError::Structure(v)=>v,_=>unreachable!("immutable native text source")})),ck::ProteinReadError::Pdb(_)=>BioPdbReadError::from_owner(owner.child(|e|match e{ck::ProteinReadError::Pdb(v)=>v,_=>unreachable!("immutable native text source")})),ck::ProteinReadError::Mmcif(_)=>BioMmcifReadError::from_owner(owner.child(|e|match e{ck::ProteinReadError::Mmcif(v)=>v,_=>unreachable!("immutable native text source")})),ck::ProteinReadError::Projection(_)=>ProteinProjectionError::from_owner(owner.child(|e|match e{ck::ProteinReadError::Projection(v)=>v,_=>unreachable!("immutable native text source")})),_=>return None}.unchecked_into()));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "JsonError | MmjsonStructureError")]
    #[derive(Clone)]
    pub type MmjsonCause;
    #[wasm_bindgen(
        typescript_type = "IoError | FromUtf8Error | BioPdbReadError | CifReadError | BioMmcifReadError | MmjsonReadError | ChemCompModelError"
    )]
    #[derive(Clone)]
    pub type BioReadCause;
}
fn mmjson_cause(owner: &ReadonlyValue<ck::MmjsonReadError>) -> Option<MmjsonCause> {
    Some(
        match owner.as_ref() {
            ck::MmjsonReadError::Parse { .. } => JsonError::from_owner(owner.child(|e| match e {
                ck::MmjsonReadError::Parse { error, .. } => error,
                _ => unreachable!("immutable JSON source"),
            })),
            ck::MmjsonReadError::Structure { .. } => {
                MmjsonStructureError::from_owner(owner.child(|e| match e {
                    ck::MmjsonReadError::Structure { error, .. } => error,
                    _ => unreachable!("immutable mmJSON structure source"),
                }))
            }
        }
        .unchecked_into(),
    )
}
enum_error!(MmjsonReadError,MmjsonReadErrorKind,ck::MmjsonReadError,{Parse=>ck::MmjsonReadError::Parse{..},Structure=>ck::MmjsonReadError::Structure{..}},{
#[wasm_bindgen(getter,js_name=sourceName)]pub fn source_name(&self)->String{match self.inner.as_ref(){ck::MmjsonReadError::Parse{source_name,..}|ck::MmjsonReadError::Structure{source_name,..}=>source_name.clone()}}
#[wasm_bindgen(getter)]pub fn error(&self)->Option<MmjsonCause>{mmjson_cause(&self.inner)}
},source=MmjsonCause,cause=mmjson_cause);
enum_error!(BioReadError,BioReadErrorKind,ck::BioReadError,{Io=>ck::BioReadError::Io{..},Utf8=>ck::BioReadError::Utf8{..},Pdb=>ck::BioReadError::Pdb(_),Cif=>ck::BioReadError::Cif(_),Mmcif=>ck::BioReadError::Mmcif(_),Mmjson=>ck::BioReadError::Mmjson(_),ChemComp=>ck::BioReadError::ChemComp(_),WrongFormat=>ck::BioReadError::WrongFormat{..},UnknownFileFormat=>ck::BioReadError::UnknownFileFormat(_)},{
#[wasm_bindgen(getter)]pub fn path(&self)->Option<String>{match self.inner.as_ref(){ck::BioReadError::Io{path,..}|ck::BioReadError::Utf8{path,..}|ck::BioReadError::UnknownFileFormat(path)=>Some(path.display().to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=sourceName)]pub fn source_name(&self)->Option<String>{match self.inner.as_ref(){ck::BioReadError::WrongFormat{source_name,..}=>Some(source_name.clone()),_=>None}}
#[wasm_bindgen(getter)]pub fn format(&self)->Option<BioCoordinateFormat>{match self.inner.as_ref(){ck::BioReadError::WrongFormat{format,..}=>Some((*format).into()),_=>None}}
},source=BioReadCause,cause=|owner:&ReadonlyValue<ck::BioReadError>|Some(match owner.as_ref(){
ck::BioReadError::Io{..}=>IoError::from_owner(owner.child(|e|match e{ck::BioReadError::Io{source,..}=>source,_=>unreachable!("immutable IO source")})),
ck::BioReadError::Utf8{..}=>FromUtf8Error::from_owner(owner.child(|e|match e{ck::BioReadError::Utf8{source,..}=>source,_=>unreachable!("immutable UTF8 source")})),
ck::BioReadError::Pdb(_)=>BioPdbReadError::from_owner(owner.child(|e|match e{ck::BioReadError::Pdb(v)=>v,_=>unreachable!("immutable format source")})),
ck::BioReadError::Cif(_)=>CifReadError::from_owner(owner.child(|e|match e{ck::BioReadError::Cif(v)=>v,_=>unreachable!("immutable format source")})),
ck::BioReadError::Mmcif(_)=>BioMmcifReadError::from_owner(owner.child(|e|match e{ck::BioReadError::Mmcif(v)=>v,_=>unreachable!("immutable format source")})),
ck::BioReadError::Mmjson(_)=>MmjsonReadError::from_owner(owner.child(|e|match e{ck::BioReadError::Mmjson(v)=>v.downcast_ref::<ck::MmjsonReadError>().expect("native Mmjson cause type"),_=>unreachable!("immutable format source")})),
ck::BioReadError::ChemComp(_)=>ChemCompModelError::from_owner(owner.child(|e|match e{ck::BioReadError::ChemComp(v)=>v.downcast_ref::<ck::ChemCompModelError>().expect("native ChemComp cause type"),_=>unreachable!("immutable format source")})),
ck::BioReadError::WrongFormat{..}|ck::BioReadError::UnknownFileFormat(_)=>return None}.unchecked_into()));
ck::__ck_configuration_BioReadParams!(crate::project_parameters, ck::BioReadParams, {}, {
    BioReadParams
});
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=fromTextWithParams)]
    pub fn from_text_with_params(
        text: Ts<input::TextInput>,
        params: BioReadParamsInput,
    ) -> Result<Self, JsValue> {
        Self::from_text(text, Some(params))
    }
    #[wasm_bindgen(js_name=fromText)]
    pub fn from_text(
        text: Ts<input::TextInput>,
        params: Option<BioReadParamsInput>,
    ) -> Result<Self, JsValue> {
        let text = input::read::<String>(text)?;
        BioReadParams::with_input(params.as_ref(), |p| {
            ck::BioStructure::from_text_with_params(&text, p)
                .map(Self)
                .map_err(BioReadError::from_rust)
        })
    }
}
#[wasm_bindgen(js_name=bioFromText)]
pub fn bio_from_text(
    text: Ts<input::TextInput>,
    params: Option<BioReadParamsInput>,
) -> Result<BioStructure, JsValue> {
    BioStructure::from_text(text, params)
}
