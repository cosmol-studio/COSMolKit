//! Structural operations retain canonical selection and native error sources.
use crate::{
    bio_hierarchy::BioStructure,
    bio_read_errors::BioStructureError,
    errors::{enum_error, error_projection, exception},
    group_values::Point3,
    input,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use std::{error::Error, sync::Arc};
use tsify::Ts;
use wasm_bindgen::{JsCast, prelude::*};
enum_error!(BioRowModelError,BioRowModelErrorKind,ck::BioRowModelError,{MissingModelNumber=>ck::BioRowModelError::MissingModelNumber},{});
enum_error!(BioRowChainError,BioRowChainErrorKind,ck::BioRowChainError,{MissingCanonicalChainName=>ck::BioRowChainError::MissingCanonicalChainName},{});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "BioRowModelError | BioRowChainError")]
    #[derive(Clone)]
    pub type BioRowTraverseErrorSource;
}
enum_error!(BioRowTraverseError,BioRowTraverseErrorKind,ck::BioRowTraverseError,{Model=>ck::BioRowTraverseError::Model(_),Chain=>ck::BioRowTraverseError::Chain(_)},{},source=BioRowTraverseErrorSource,cause=|owner:&ReadonlyValue<ck::BioRowTraverseError>|Some(match owner.as_ref(){ck::BioRowTraverseError::Model(_)=>BioRowModelError::from_owner(owner.child(|e|match e{ck::BioRowTraverseError::Model(v)=>v,_=>unreachable!("immutable structural error source")})),ck::BioRowTraverseError::Chain(_)=>BioRowChainError::from_owner(owner.child(|e|match e{ck::BioRowTraverseError::Chain(v)=>v,_=>unreachable!("immutable structural error source")}))}.unchecked_into()));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "BioStructureError | BioRowTraverseError")]
    #[derive(Clone)]
    pub type BioSelectionCopyCauseSource;
}
enum_error!(BioSelectionCopyCause,BioSelectionCopyCauseKind,ck::BioSelectionCopyCause,{Structure=>ck::BioSelectionCopyCause::Structure(_),Traverse=>ck::BioSelectionCopyCause::Traverse(_)},{},source=BioSelectionCopyCauseSource,cause=|owner:&ReadonlyValue<ck::BioSelectionCopyCause>|Some(match owner.as_ref(){ck::BioSelectionCopyCause::Structure(_)=>BioStructureError::from_owner(owner.child(|e|match e{ck::BioSelectionCopyCause::Structure(v)=>v,_=>unreachable!("immutable structural error source")})),ck::BioSelectionCopyCause::Traverse(_)=>BioRowTraverseError::from_owner(owner.child(|e|match e{ck::BioSelectionCopyCause::Traverse(v)=>v,_=>unreachable!("immutable structural error source")}))}.unchecked_into()));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "BioStructureError | ProteinProjectionError | BioSelectionCopyError"
    )]
    #[derive(Clone)]
    pub type BioOperationErrorSource;
}
enum_error!(BioOperationError,BioOperationErrorKind,ck::BioOperationError,{Structure=>ck::BioOperationError::Structure(_),Protein=>ck::BioOperationError::Protein(_),Selection=>ck::BioOperationError::Selection(_)},{},source=BioOperationErrorSource,cause=|owner:&ReadonlyValue<ck::BioOperationError>|Some(match owner.as_ref(){ck::BioOperationError::Structure(_)=>BioStructureError::from_owner(owner.child(|e|match e{ck::BioOperationError::Structure(v)=>v,_=>unreachable!("immutable structural error source")})),ck::BioOperationError::Protein(_)=>ProteinProjectionError::from_owner(owner.child(|e|match e{ck::BioOperationError::Protein(v)=>v,_=>unreachable!("immutable structural error source")})),ck::BioOperationError::Selection(_)=>BioSelectionCopyError::from_owner(owner.child(|e|match e{ck::BioOperationError::Selection(v)=>v,_=>unreachable!("immutable structural error source")}))}.unchecked_into()));
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "BioStructureError")]
    #[derive(Clone)]
    pub type ProteinProjectionCause;
    #[wasm_bindgen(typescript_type = "BioRowTraverseError")]
    #[derive(Clone)]
    pub type BioSelectionMatchCause;
    #[wasm_bindgen(typescript_type = "BioSelectionCopyCause")]
    pub type BioSelectionCopyPayload;
}
enum_error!(ProteinProjectionError,ProteinProjectionErrorKind,ck::ProteinProjectionError,{Structure=>ck::ProteinProjectionError::Structure(_),MissingEntityMapping=>ck::ProteinProjectionError::MissingEntityMapping{..},NonAminoAcidResidue=>ck::ProteinProjectionError::NonAminoAcidResidue{..}},{
#[wasm_bindgen(getter)]pub fn source(&self)->Option<u32>{match self.inner.as_ref(){ck::ProteinProjectionError::MissingEntityMapping{source}=>Some(source.value()),_=>None}}
#[wasm_bindgen(getter)]pub fn index(&self)->Option<usize>{match self.inner.as_ref(){ck::ProteinProjectionError::NonAminoAcidResidue{index}=>Some(*index),_=>None}}
},source=ProteinProjectionCause,cause=|owner:&ReadonlyValue<ck::ProteinProjectionError>|match owner.as_ref(){ck::ProteinProjectionError::Structure(_)=>Some(BioStructureError::from_owner(owner.child(|e|match e{ck::ProteinProjectionError::Structure(v)=>v,_=>unreachable!("immutable protein structure error")})).unchecked_into()),_=>None});
error_projection!(
    BioSelectionMatchError,
    ck::BioSelectionMatchError,
    {},
    source = BioSelectionMatchCause,
    cause = |owner: &ReadonlyValue<ck::BioSelectionMatchError>| Some(
        BioRowTraverseError::from_owner(owner.child(|e| {
            Error::source(e)
                .expect("native selection source")
                .downcast_ref::<ck::BioRowTraverseError>()
                .expect("native traversal source")
        }))
        .unchecked_into()
    )
);
error_projection!(
    BioSelectionCopyError,
    ck::BioSelectionCopyError,
    {
        #[wasm_bindgen(getter)]
        pub fn cause(&self) -> BioSelectionCopyPayload {
            BioSelectionCopyCause::from_owner(self.inner.child(ck::BioSelectionCopyError::cause))
                .unchecked_into()
        }
    },
    source = BioSelectionCopyCauseSource,
    cause = |owner: &ReadonlyValue<ck::BioSelectionCopyError>| Some(
        match owner.cause() {
            ck::BioSelectionCopyCause::Structure(_) =>
                BioStructureError::from_owner(owner.child(|e| match e.cause() {
                    ck::BioSelectionCopyCause::Structure(v) => v,
                    _ => unreachable!("immutable copy structure error"),
                })),
            ck::BioSelectionCopyCause::Traverse(_) =>
                BioRowTraverseError::from_owner(owner.child(|e| match e.cause() {
                    ck::BioSelectionCopyCause::Traverse(v) => v,
                    _ => unreachable!("immutable copy traversal error"),
                })),
        }
        .unchecked_into()
    )
);
enum_error!(BioSelectionParseError,BioSelectionParseErrorKind,ck::BioSelectionParseError,{Syntax=>ck::BioSelectionParseError::Syntax(_),SeqidRange=>ck::BioSelectionParseError::SeqidRange(_)},{
#[wasm_bindgen(getter)]pub fn cid(&self)->Option<String>{match self.inner.as_ref(){ck::BioSelectionParseError::Syntax(v)=>Some(v.cid().into()),_=>None}}
#[wasm_bindgen(getter)]pub fn pos(&self)->Option<usize>{match self.inner.as_ref(){ck::BioSelectionParseError::Syntax(v)=>Some(v.pos()),_=>None}}
#[wasm_bindgen(getter)]pub fn info(&self)->Option<String>{match self.inner.as_ref(){ck::BioSelectionParseError::Syntax(v)=>v.info().map(str::to_owned),_=>None}}
});
#[wasm_bindgen]
pub struct BioSelection(Arc<ck::BioSelection>);
crate::fields::value_field!(
    Arc<ck::BioSelection>,
    BioSelection,
    BioSelectionInput,
    "BioSelection"
);
#[wasm_bindgen]
impl BioSelection {
    #[wasm_bindgen(js_name=fromCid)]
    pub fn from_cid(cid: Ts<input::TextInput>) -> Result<Self, JsValue> {
        ck::BioSelection::from_cid(&input::read::<String>(cid)?)
            .map(|v| Self(Arc::new(v)))
            .map_err(BioSelectionParseError::from_rust)
    }
    #[wasm_bindgen(js_name=toCid)]
    pub fn to_cid(&self) -> String {
        self.0.to_cid()
    }
}
#[wasm_bindgen]
impl BioStructure {
    #[wasm_bindgen(js_name=selectedAtomIds)]
    pub fn selected_atom_ids(&self, selection: BioSelectionInput) -> Result<Vec<u32>, JsValue> {
        BioSelection::with_ref(&selection, |s| {
            self.0
                .selected_atom_ids(s)
                .map(|v| v.into_iter().map(|id| id.value()).collect())
                .map_err(BioSelectionMatchError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withSelection)]
    pub fn with_selection(&self, selection: BioSelectionInput) -> Result<Self, JsValue> {
        BioSelection::with_ref(&selection, |s| {
            self.0
                .with_selection(s)
                .map(Self)
                .map_err(BioOperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=retainSelection_)]
    pub fn retain_selection_(&mut self, selection: BioSelectionInput) -> Result<(), JsValue> {
        BioSelection::with_ref(&selection, |s| {
            self.0
                .retain_selection_(s)
                .map_err(BioOperationError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=withTranslatedCoordinates)]
    pub fn with_translated_coordinates(&self, offset: Ts<Point3>) -> Result<Self, JsValue> {
        self.0
            .with_translated_coordinates(input::decode(offset)?.0)
            .map(Self)
            .map_err(BioOperationError::from_rust)
    }
    #[wasm_bindgen(js_name=translate_)]
    pub fn translate_(&mut self, offset: Ts<Point3>) -> Result<(), JsValue> {
        self.0
            .translate_(input::decode(offset)?.0)
            .map_err(BioOperationError::from_rust)
    }
}
