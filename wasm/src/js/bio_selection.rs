use crate::bio_metadata::BioTransform;
// Thin BIO selection and COW transformation projections.
use crate::{
    alignment_values::{set, source_error},
    bio_readers::BioStructure,
    host_values::type_error,
};
use cosmolkit_wasm::rust as ck;
use std::{error::Error as RustError, rc::Rc};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BioSelection {
    pub(crate) inner: ck::BioSelection,
}
#[wasm_bindgen]
impl BioSelection {
    #[wasm_bindgen(js_name=fromCid)]
    pub fn from_cid(cid: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::BioSelection::from_cid(cid).map(|inner| Self { inner })
        ck::BioSelection::from_cid(cid)
            .map(|inner| Self { inner })
            .map_err(|e| parse_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toCid)]
    pub fn to_cid(&self) -> String {
        self.inner.to_cid()
    }
}
pub(crate) use crate::host_values::point;
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=selectedAtomIds)]
    pub fn selected_atom_ids(&self, selection: &BioSelection) -> Result<Vec<u32>, JsValue> {
        self.inner
            .selected_atom_ids(&selection.inner)
            .map(|ids| ids.into_iter().map(|id| id.value()).collect())
            .map_err(|e| match_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withSelection)]
    pub fn with_selection(&self, selection: &BioSelection) -> Result<Self, JsValue> {
        self.inner
            .with_selection(&selection.inner)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=retainSelection_)]
    pub fn retain_selection_(&mut self, selection: &BioSelection) -> Result<(), JsValue> {
        // COSMolKit❗✔️: Arc::make_mut(&mut self.inner).retain_selection_(&selection.inner)
        // The facade owns COW and atomic commit. Existing row snapshots retain their owner.
        Rc::make_mut(&mut self.inner)
            .retain_selection_(&selection.inner)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withTranslatedCoordinates)]
    pub fn with_translated_coordinates(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array")] offset: JsValue,
    ) -> Result<Self, JsValue> {
        self.inner
            .with_translated_coordinates(point(&offset)?)
            .map(|inner| Self {
                inner: Rc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=translate_)]
    pub fn translate_(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Float64Array")] offset: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: Arc::make_mut(&mut self.inner).translate_(offset)
        // The validated public operation is the sole coordinate write boundary.
        let offset = point(&offset)?;
        Rc::make_mut(&mut self.inner)
            .translate_(offset)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl BioTransform {
    pub fn approx(
        &self,
        other: &crate::bio_metadata::BioTransform,
        #[wasm_bindgen(unchecked_param_type = "number")] epsilon: JsValue,
    ) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.approx(&other.inner, epsilon)
        let epsilon = epsilon.as_f64().ok_or_else(|| type_error("epsilon"))?;
        Ok(self.inner.approx(&other.inner, epsilon))
    }
}

#[wasm_bindgen]
pub struct BioSelectionParseError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioSelectionParseError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string")]
    pub fn variant(&self) -> JsValue {
        field(&self.context, "variant")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn cid(&self) -> JsValue {
        field(&self.context, "cid")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn pos(&self) -> JsValue {
        field(&self.context, "pos")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn info(&self) -> JsValue {
        field(&self.context, "info")
    }
}
fn make_bioselectionparseerror(
    source: &ck::BioSelectionParseError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioSelectionParseError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        js_sys::Reflect::set(&e, &key, &js_sys::Reflect::get(&context, &key)?)?;
    }
    set(
        &e,
        "detail",
        BioSelectionParseError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
#[wasm_bindgen]
pub struct BioSelectionMatchError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioSelectionMatchError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
fn make_bioselectionmatcherror(
    source: &ck::BioSelectionMatchError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioSelectionMatchError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        js_sys::Reflect::set(&e, &key, &js_sys::Reflect::get(&context, &key)?)?;
    }
    set(
        &e,
        "detail",
        BioSelectionMatchError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
#[wasm_bindgen]
pub struct BioOperationError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioOperationError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
fn make_biooperationerror(
    source: &ck::BioOperationError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioOperationError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        js_sys::Reflect::set(&e, &key, &js_sys::Reflect::get(&context, &key)?)?;
    }
    set(
        &e,
        "detail",
        BioOperationError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
#[wasm_bindgen]
pub struct ProteinProjectionError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl ProteinProjectionError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn index(&self) -> JsValue {
        field(&self.context, "index")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn source(&self) -> JsValue {
        field(&self.context, "source")
    }
}
fn make_proteinprojectionerror(
    source: &ck::ProteinProjectionError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("ProteinProjectionError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        js_sys::Reflect::set(&e, &key, &js_sys::Reflect::get(&context, &key)?)?;
    }
    set(
        &e,
        "detail",
        ProteinProjectionError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
#[wasm_bindgen]
pub struct BioSelectionCopyCause {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioSelectionCopyCause {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
fn make_bioselectioncopycause(
    source: &ck::BioSelectionCopyCause,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioSelectionCopyCause");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        js_sys::Reflect::set(&e, &key, &js_sys::Reflect::get(&context, &key)?)?;
    }
    set(
        &e,
        "detail",
        BioSelectionCopyCause {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
fn field(context: &JsValue, key: &str) -> JsValue {
    let v = js_sys::Reflect::get(context, &key.into()).unwrap_or(JsValue::NULL);
    if v.is_undefined() { JsValue::NULL } else { v }
}
pub(crate) fn parse_error(source: &ck::BioSelectionParseError) -> Result<JsValue, JsValue> {
    let context: JsValue = js_sys::Object::new().into();
    match source {
        ck::BioSelectionParseError::Syntax(s) => {
            set(&context, "variant", "Syntax".into())?;
            set(&context, "cid", s.cid().into())?;
            set(&context, "pos", (s.pos() as f64).into())?;
            set(
                &context,
                "info",
                s.info().map_or(JsValue::NULL, JsValue::from),
            )?;
        }
        ck::BioSelectionParseError::SeqidRange(_) => {
            set(&context, "variant", "SeqidRange".into())?;
        }
    }
    make_bioselectionparseerror(source, "Parse", context)
}
pub(crate) fn match_error(source: &ck::BioSelectionMatchError) -> Result<JsValue, JsValue> {
    make_bioselectionmatcherror(source, "SelectionMatch", js_sys::Object::new().into())
}
pub(crate) fn operation_error(source: &ck::BioOperationError) -> Result<JsValue, JsValue> {
    make_biooperationerror(
        source,
        match source {
            ck::BioOperationError::Structure(_) => "Structure",
            ck::BioOperationError::Protein(_) => "Protein",
            ck::BioOperationError::Selection(_) => "Selection",
        },
        js_sys::Object::new().into(),
    )
}
pub(crate) fn protein_projection_error(
    source: &ck::ProteinProjectionError,
) -> Result<JsValue, JsValue> {
    let context: JsValue = js_sys::Object::new().into();
    let kind = match source {
        ck::ProteinProjectionError::Structure(_) => "Structure",
        ck::ProteinProjectionError::MissingEntityMapping { source } => {
            set(&context, "source", source.value().into())?;
            "MissingEntityMapping"
        }
        ck::ProteinProjectionError::NonAminoAcidResidue { index } => {
            set(&context, "index", (*index as f64).into())?;
            "NonAminoAcidResidue"
        }
    };
    make_proteinprojectionerror(source, kind, context)
}
pub(crate) fn copy_cause_error(source: &ck::BioSelectionCopyCause) -> Result<JsValue, JsValue> {
    make_bioselectioncopycause(
        source,
        match source {
            ck::BioSelectionCopyCause::Structure(_) => "Structure",
            ck::BioSelectionCopyCause::Traverse(_) => "Traverse",
        },
        js_sys::Object::new().into(),
    )
}
#[wasm_bindgen]
pub struct BioSelectionCopyError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl BioSelectionCopyError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(unchecked_return_type = "BioSelectionCopyCause")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
pub(crate) fn copy_error(source: &ck::BioSelectionCopyError) -> Result<JsValue, JsValue> {
    let kind = match source.cause() {
        ck::BioSelectionCopyCause::Structure(_) => "Structure",
        ck::BioSelectionCopyCause::Traverse(_) => "Traverse",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let typed_cause = BioSelectionCopyCause {
        kind: kind.into(),
        message: source.cause().to_string(),
        cause: cause.clone(),
        context: js_sys::Object::new().into(),
    };
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioSelectionCopyError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    set(
        &e,
        "detail",
        BioSelectionCopyError {
            kind: kind.into(),
            message: source.to_string(),
            cause: typed_cause.into(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
