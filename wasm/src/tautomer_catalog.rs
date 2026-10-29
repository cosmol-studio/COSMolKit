//! Catalog construction remains in the native tautomer owner.
use crate::{
    errors::enum_error, fields::Field, input, storage::NativeField, storage::ReadonlyValue,
    tautomer_calls::TautomerParams,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "SmartsParseError | QueryCompileError")]
    #[derive(Clone)]
    pub type TransformCause;
    #[wasm_bindgen(typescript_type = "IoError | FromUtf8Error | TautomerTransformError")]
    #[derive(Clone)]
    pub type CatalogCause;
}
fn transform_cause(owner: &ReadonlyValue<ck::TautomerTransformError>) -> Option<TransformCause> {
    Some(match owner.as_ref() {
        ck::TautomerTransformError::ChargeSymbolNotRecognised => return None,
        ck::TautomerTransformError::CannotParseSmarts { .. } => {
            crate::smarts_parse_errors::SmartsParseError::from_owner(owner.child(|e| match e {
                ck::TautomerTransformError::CannotParseSmarts { source, .. } => source,
                _ => unreachable!(),
            }))
            .unchecked_into()
        }
        ck::TautomerTransformError::QueryCompile(_) => {
            crate::query_compile_error::QueryCompileError::from_owner(owner.child(|e| match e {
                ck::TautomerTransformError::QueryCompile(source) => source,
                _ => unreachable!(),
            }))
            .unchecked_into()
        }
    })
}
enum_error!(TautomerTransformError,TautomerTransformErrorKind,ck::TautomerTransformError,{
    ChargeSymbolNotRecognised=>ck::TautomerTransformError::ChargeSymbolNotRecognised,CannotParseSmarts=>ck::TautomerTransformError::CannotParseSmarts{..},QueryCompile=>ck::TautomerTransformError::QueryCompile(_)
},{
    #[wasm_bindgen(getter)]pub fn smarts(&self)->Option<String>{match self.inner.as_ref(){ck::TautomerTransformError::CannotParseSmarts{smarts,..}=>Some(smarts.clone()),_=>None}}
},source=TransformCause,cause=transform_cause);
fn catalog_cause(owner: &ReadonlyValue<ck::TautomerCatalogError>) -> Option<CatalogCause> {
    Some(match owner.as_ref() {
        ck::TautomerCatalogError::BadInputFile { .. }
        | ck::TautomerCatalogError::BadStreamContents(_) => {
            crate::external_errors::IoError::from_owner(owner.child(|e| match e {
                ck::TautomerCatalogError::BadInputFile { source, .. }
                | ck::TautomerCatalogError::BadStreamContents(source) => source,
                _ => unreachable!(),
            }))
            .unchecked_into()
        }
        ck::TautomerCatalogError::InvalidUtf8(_) => {
            crate::external_errors::FromUtf8Error::from_owner(owner.child(|e| match e {
                ck::TautomerCatalogError::InvalidUtf8(source) => source,
                _ => unreachable!(),
            }))
            .unchecked_into()
        }
        ck::TautomerCatalogError::Transform(_) => {
            TautomerTransformError::from_owner(owner.child(|e| match e {
                ck::TautomerCatalogError::Transform(source) => source,
                _ => unreachable!(),
            }))
            .unchecked_into()
        }
        _ => return None,
    })
}
enum_error!(TautomerCatalogError,TautomerCatalogErrorKind,ck::TautomerCatalogError,{
    BadInputFile=>ck::TautomerCatalogError::BadInputFile{..},BadStreamContents=>ck::TautomerCatalogError::BadStreamContents(_),InvalidUtf8=>ck::TautomerCatalogError::InvalidUtf8(_),Transform=>ck::TautomerCatalogError::Transform(_),TransformIndexOutOfRange=>ck::TautomerCatalogError::TransformIndexOutOfRange{..},DeserializationUnderConstruction=>ck::TautomerCatalogError::DeserializationUnderConstruction
},{
    #[wasm_bindgen(getter)]pub fn path(&self)->Option<String>{match self.inner.as_ref(){ck::TautomerCatalogError::BadInputFile{path,..}=>Some(path.display().to_string()),_=>None}}
    #[wasm_bindgen(getter)]pub fn index(&self)->Option<usize>{match self.inner.as_ref(){ck::TautomerCatalogError::TransformIndexOutOfRange{index,..}=>Some(*index),_=>None}}
    #[wasm_bindgen(getter)]pub fn len(&self)->Option<usize>{match self.inner.as_ref(){ck::TautomerCatalogError::TransformIndexOutOfRange{len,..}=>Some(*len),_=>None}}
},source=CatalogCause,cause=catalog_cause);

#[derive(serde::Deserialize, tsify::Tsify)]
#[serde(transparent)]
pub struct TransformRows(
    #[tsify(type = "ReadonlyArray<readonly [string, string, string, string]>")] Vec<Vec<String>>,
);
#[wasm_bindgen]
impl TautomerParams {
    pub fn v1() -> Result<Self, JsValue> {
        let value = ck::TautomerParams::v1().map_err(TautomerCatalogError::from_rust)?;
        ck::TautomerParams::get(NativeField::owned(value))
    }
    #[wasm_bindgen(js_name=fromTransformData)]
    pub fn from_transform_data(data: tsify::Ts<TransformRows>) -> Result<Self, JsValue> {
        let rows = input::decode(data)?
            .0
            .into_iter()
            .map(|row| {
                input::fixed_array::<String, 4, serde::de::value::Error>(row)
                    .map_err(|e| JsValue::from(js_sys::TypeError::new(&e.to_string())))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let borrowed = rows
            .iter()
            .map(|r| (r[0].as_str(), r[1].as_str(), r[2].as_str(), r[3].as_str()))
            .collect::<Vec<_>>();
        let value = ck::TautomerParams::from_transform_data(&borrowed)
            .map_err(TautomerCatalogError::from_rust)?;
        ck::TautomerParams::get(NativeField::owned(value))
    }
}
