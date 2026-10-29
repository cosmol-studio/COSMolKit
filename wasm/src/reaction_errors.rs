//! Exact native reaction parsing, template and validation failures.
use crate::external_errors::Utf8Error;
use crate::model_errors::AtomPropertyError;
use crate::path_errors::PathError;
use crate::property_errors::PropertyIntReadError;
use crate::query_model_errors::QueryGraphError;
use crate::smarts_errors::SmartsWriteError;
use crate::smarts_parse_errors::{CxQueryLoweringError, SmartsParseError};
use crate::smiles_parse_errors::SmilesParseError;
use crate::stereo_order_errors::StereoOrderError;
use crate::{
    errors::{enum_error, error_projection, exception},
    reaction_values::*,
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
error_projection!(CxParseError, ck::CxParseError, {
    #[wasm_bindgen(getter)]
    pub fn offset(&self) -> usize {
        self.inner.offset
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.message.clone()
    }
});
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "QueryGraphError")]
    #[derive(Clone)]
    pub type ReactionModelErrorCause;
}
fn reactionmodelerror_cause(
    owner: &ReadonlyValue<ck::ReactionModelError>,
) -> Option<ReactionModelErrorCause> {
    match owner.as_ref() {
        ck::ReactionModelError::Template { .. } => Some(
            QueryGraphError::from_owner(owner.child(|e| match e {
                ck::ReactionModelError::Template { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ReactionModelError,ReactionModelErrorKind,ck::ReactionModelError,{Template=>ck::ReactionModelError::Template{..},TemplateIndex=>ck::ReactionModelError::TemplateIndex{..}},{
#[wasm_bindgen(getter,js_name=role)]pub fn role(&self)->Option<ReactionRole>{match self.inner.as_ref(){ck::ReactionModelError::Template{role,..}=>Some((*role).into()),ck::ReactionModelError::TemplateIndex{role,..}=>Some((*role).into()),_=>None}}
#[wasm_bindgen(getter,js_name=template)]pub fn template(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionModelError::Template{template,..}=>Some(*template as usize),_=>None}}
#[wasm_bindgen(getter,js_name=index)]pub fn index(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionModelError::TemplateIndex{index,..}=>Some(*index as usize),_=>None}}
#[wasm_bindgen(getter,js_name=count)]pub fn count(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionModelError::TemplateIndex{count,..}=>Some(*count as usize),_=>None}}
},source=ReactionModelErrorCause,cause=reactionmodelerror_cause);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "PropertyIntReadError | AtomPropertyError")]
    #[derive(Clone)]
    pub type ReactionValidationErrorCause;
}
fn reactionvalidationerror_cause(
    owner: &ReadonlyValue<ck::ReactionValidationError>,
) -> Option<ReactionValidationErrorCause> {
    match owner.as_ref() {
        ck::ReactionValidationError::Property { .. } => Some(
            PropertyIntReadError::from_owner(owner.child(|e| match e {
                ck::ReactionValidationError::Property { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionValidationError::Annotation { .. } => Some(
            AtomPropertyError::from_owner(owner.child(|e| match e {
                ck::ReactionValidationError::Annotation { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ReactionValidationError,ReactionValidationErrorKind,ck::ReactionValidationError,{Property=>ck::ReactionValidationError::Property{..},MapOverflow=>ck::ReactionValidationError::MapOverflow{..},AtomBounds=>ck::ReactionValidationError::AtomBounds{..},MissingReactingAtom=>ck::ReactionValidationError::MissingReactingAtom{..},QueryArithmetic=>ck::ReactionValidationError::QueryArithmetic{..},Annotation=>ck::ReactionValidationError::Annotation{..}},{
#[wasm_bindgen(getter,js_name=role)]pub fn role(&self)->Option<ReactionRole>{match self.inner.as_ref(){ck::ReactionValidationError::Property{role,..}=>Some((*role).into()),ck::ReactionValidationError::MapOverflow{role,..}=>Some((*role).into()),ck::ReactionValidationError::AtomBounds{role,..}=>Some((*role).into()),_=>None}}
#[wasm_bindgen(getter,js_name=template)]pub fn template(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionValidationError::Property{template,..}=>Some(*template as usize),ck::ReactionValidationError::MapOverflow{template,..}=>Some(*template as usize),ck::ReactionValidationError::AtomBounds{template,..}=>Some(*template as usize),ck::ReactionValidationError::MissingReactingAtom{template,..}=>Some(*template as usize),ck::ReactionValidationError::QueryArithmetic{template,..}=>Some(*template as usize),ck::ReactionValidationError::Annotation{template,..}=>Some(*template as usize),_=>None}}
#[wasm_bindgen(getter,js_name=atom)]pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionValidationError::Property{atom,..}=>Some(atom.index()),ck::ReactionValidationError::MapOverflow{atom,..}=>Some(atom.index()),ck::ReactionValidationError::AtomBounds{atom,..}=>Some(atom.index()),ck::ReactionValidationError::MissingReactingAtom{atom,..}=>Some(atom.index()),ck::ReactionValidationError::QueryArithmetic{atom,..}=>Some(atom.index()),ck::ReactionValidationError::Annotation{atom,..}=>Some(atom.index()),_=>None}}
#[wasm_bindgen(getter,js_name=property)]pub fn property(&self)->Option<String>{match self.inner.as_ref(){ck::ReactionValidationError::Property{property,..}=>Some(property.to_string()),ck::ReactionValidationError::MapOverflow{property,..}=>Some(property.to_string()),ck::ReactionValidationError::Annotation{property,..}=>Some(property.to_string()),_=>None}}
#[wasm_bindgen(getter,js_name=value)]pub fn value(&self)->Option<u32>{match self.inner.as_ref(){ck::ReactionValidationError::MapOverflow{value,..}=>Some(*value),_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionValidationError::AtomBounds{atom_count,..}=>Some(*atom_count as usize),_=>None}}
#[wasm_bindgen(getter,js_name=map)]pub fn map(&self)->Option<i32>{match self.inner.as_ref(){ck::ReactionValidationError::MissingReactingAtom{map,..}=>Some(*map),_=>None}}
},source=ReactionValidationErrorCause,cause=reactionvalidationerror_cause);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReactionValidationError")]
    #[derive(Clone)]
    pub type ReactionInitializationErrorCause;
}
fn reactioninitializationerror_cause(
    owner: &ReadonlyValue<ck::ReactionInitializationError>,
) -> Option<ReactionInitializationErrorCause> {
    match owner.as_ref() {
        ck::ReactionInitializationError::Validation { .. } => Some(
            ReactionValidationError::from_owner(owner.child(|e| match e {
                ck::ReactionInitializationError::Validation { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ReactionInitializationError,ReactionInitializationErrorKind,ck::ReactionInitializationError,{Validation=>ck::ReactionInitializationError::Validation{..},Invalid=>ck::ReactionInitializationError::Invalid{..}},{
#[wasm_bindgen(getter,js_name=report)]pub fn report(&self)->Option<ReactionValidationReport>{match self.inner.as_ref(){ck::ReactionInitializationError::Invalid{report,..}=>Some(ReactionValidationReport(report.clone())),_=>None}}
},source=ReactionInitializationErrorCause,cause=reactioninitializationerror_cause);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "SmartsParseError | ReactionValidationError | SmilesParseError | Utf8Error | QueryGraphError | CxParseError | CxQueryLoweringError | StereoOrderError"
    )]
    #[derive(Clone)]
    pub type ReactionParseErrorCause;
}
fn reactionparseerror_cause(
    owner: &ReadonlyValue<ck::ReactionParseError>,
) -> Option<ReactionParseErrorCause> {
    match owner.as_ref() {
        ck::ReactionParseError::ParserCleanup { .. } => Some(
            SmartsParseError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::ParserCleanup { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::TemplateProperty { .. } => Some(
            ReactionValidationError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::TemplateProperty { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::Smarts { .. } => Some(
            SmartsParseError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::Smarts { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::Smiles { .. } => Some(
            SmilesParseError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::Smiles { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::ComponentEncoding { .. } => Some(
            Utf8Error::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::ComponentEncoding { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::Model { .. } => Some(
            QueryGraphError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::Model { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::AgentFragments { .. } => Some(
            SmartsParseError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::AgentFragments { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::CxParse { .. } => Some(
            CxParseError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::CxParse { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::CxLowering { .. } => Some(
            CxQueryLoweringError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::CxLowering { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionParseError::StereoOrder { .. } => Some(
            StereoOrderError::from_owner(owner.child(|e| match e {
                ck::ReactionParseError::StereoOrder { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ReactionParseError,ReactionParseErrorKind,ck::ReactionParseError,{ParserCleanup=>ck::ReactionParseError::ParserCleanup{..},TemplateAtomBounds=>ck::ReactionParseError::TemplateAtomBounds{..},StereoDegreeArithmetic=>ck::ReactionParseError::StereoDegreeArithmetic{..},StereoInversionFlag=>ck::ReactionParseError::StereoInversionFlag{..},TemplateProperty=>ck::ReactionParseError::TemplateProperty{..},Separators=>ck::ReactionParseError::Separators{..},MultiStep=>ck::ReactionParseError::MultiStep{..},ComponentBounds=>ck::ReactionParseError::ComponentBounds{..},OffsetOverflow=>ck::ReactionParseError::OffsetOverflow,Smarts=>ck::ReactionParseError::Smarts{..},Smiles=>ck::ReactionParseError::Smiles{..},ComponentEncoding=>ck::ReactionParseError::ComponentEncoding{..},Model=>ck::ReactionParseError::Model{..},AgentFragments=>ck::ReactionParseError::AgentFragments{..},CxParse=>ck::ReactionParseError::CxParse{..},CxLowering=>ck::ReactionParseError::CxLowering{..},StereoOrder=>ck::ReactionParseError::StereoOrder{..}},{
#[wasm_bindgen(getter,js_name=role)]pub fn role(&self)->Option<ReactionRole>{match self.inner.as_ref(){ck::ReactionParseError::ParserCleanup{role,..}=>Some((*role).into()),ck::ReactionParseError::TemplateAtomBounds{role,..}=>Some((*role).into()),ck::ReactionParseError::Smarts{role,..}=>Some((*role).into()),ck::ReactionParseError::Smiles{role,..}=>Some((*role).into()),ck::ReactionParseError::ComponentEncoding{role,..}=>Some((*role).into()),ck::ReactionParseError::Model{role,..}=>Some((*role).into()),ck::ReactionParseError::CxParse{role,..}=>Some((*role).into()),ck::ReactionParseError::CxLowering{role,..}=>Some((*role).into()),_=>None}}
#[wasm_bindgen(getter,js_name=template)]pub fn template(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::ParserCleanup{template,..}=>Some(*template as usize),ck::ReactionParseError::TemplateAtomBounds{template,..}=>Some(*template as usize),ck::ReactionParseError::Smarts{template,..}=>Some(*template as usize),ck::ReactionParseError::Smiles{template,..}=>Some(*template as usize),ck::ReactionParseError::ComponentEncoding{template,..}=>Some(*template as usize),ck::ReactionParseError::Model{template,..}=>Some(*template as usize),ck::ReactionParseError::CxParse{template,..}=>Some(*template as usize),ck::ReactionParseError::CxLowering{template,..}=>Some(*template as usize),_=>None}}
#[wasm_bindgen(getter,js_name=atom)]pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::TemplateAtomBounds{atom,..}=>Some(atom.index()),ck::ReactionParseError::StereoInversionFlag{atom,..}=>Some(atom.index()),ck::ReactionParseError::StereoOrder{atom,..}=>Some(atom.index()),_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)]pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::TemplateAtomBounds{atom_count,..}=>Some(*atom_count as usize),_=>None}}
#[wasm_bindgen(getter,js_name=reactant)]pub fn reactant(&self)->Option<u32>{match self.inner.as_ref(){ck::ReactionParseError::StereoDegreeArithmetic{reactant,..}=>Some(*reactant),_=>None}}
#[wasm_bindgen(getter,js_name=product)]pub fn product(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::StereoDegreeArithmetic{product,..}=>Some(*product as usize),ck::ReactionParseError::StereoInversionFlag{product,..}=>Some(*product as usize),ck::ReactionParseError::StereoOrder{product,..}=>Some(*product as usize),_=>None}}
#[wasm_bindgen(getter,js_name=flag)]pub fn flag(&self)->Option<i32>{match self.inner.as_ref(){ck::ReactionParseError::StereoInversionFlag{flag,..}=>*flag,_=>None}}
#[wasm_bindgen(getter,js_name=count)]pub fn count(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::Separators{count,..}=>Some(*count as usize),ck::ReactionParseError::MultiStep{count,..}=>Some(*count as usize),_=>None}}
#[wasm_bindgen(getter,js_name=start)]pub fn start(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::ComponentBounds{start,..}=>Some(*start as usize),_=>None}}
#[wasm_bindgen(getter,js_name=end)]pub fn end(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::ComponentBounds{end,..}=>Some(*end as usize),_=>None}}
#[wasm_bindgen(getter,js_name=text)]pub fn text(&self)->Result<Option<String>,JsValue>{Ok(match self.inner.as_ref(){ck::ReactionParseError::Smarts{text,..}=>Some(crate::property_errors::text(text)?),ck::ReactionParseError::Smiles{text,..}=>Some(crate::property_errors::text(text)?),ck::ReactionParseError::ComponentEncoding{text,..}=>Some(crate::property_errors::text(text)?),_=>None})}
#[wasm_bindgen(getter,js_name=startAtom)]pub fn start_atom(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::CxParse{start_atom,..}=>Some(*start_atom as usize),ck::ReactionParseError::CxLowering{start_atom,..}=>Some(*start_atom as usize),_=>None}}
#[wasm_bindgen(getter,js_name=startBond)]pub fn start_bond(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionParseError::CxParse{start_bond,..}=>Some(*start_bond as usize),ck::ReactionParseError::CxLowering{start_bond,..}=>Some(*start_bond as usize),_=>None}}
},source=ReactionParseErrorCause,cause=reactionparseerror_cause);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "SmartsWriteError | PathError")]
    #[derive(Clone)]
    pub type ReactionWriteErrorCause;
}
fn reactionwriteerror_cause(
    owner: &ReadonlyValue<ck::ReactionWriteError>,
) -> Option<ReactionWriteErrorCause> {
    match owner.as_ref() {
        ck::ReactionWriteError::Cx(_) => Some(
            SmartsWriteError::from_owner(owner.child(|e| match e {
                ck::ReactionWriteError::Cx(source) => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionWriteError::Template { .. } => Some(
            SmartsWriteError::from_owner(owner.child(|e| match e {
                ck::ReactionWriteError::Template { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        ck::ReactionWriteError::Connectivity { .. } => Some(
            PathError::from_owner(owner.child(|e| match e {
                ck::ReactionWriteError::Connectivity { source, .. } => source,
                _ => unreachable!("immutable reaction cause"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ReactionWriteError,ReactionWriteErrorKind,ck::ReactionWriteError,{Cx=>ck::ReactionWriteError::Cx(_),Template=>ck::ReactionWriteError::Template{..},Connectivity=>ck::ReactionWriteError::Connectivity{..},CoordinateSelectionArity=>ck::ReactionWriteError::CoordinateSelectionArity{..}},{
#[wasm_bindgen(getter,js_name=role)]pub fn role(&self)->Option<ReactionRole>{match self.inner.as_ref(){ck::ReactionWriteError::Template{role,..}=>Some((*role).into()),ck::ReactionWriteError::Connectivity{role,..}=>Some((*role).into()),_=>None}}
#[wasm_bindgen(getter,js_name=template)]pub fn template(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionWriteError::Template{template,..}=>Some(*template as usize),ck::ReactionWriteError::Connectivity{template,..}=>Some(*template as usize),_=>None}}
#[wasm_bindgen(getter,js_name=expected)]pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionWriteError::CoordinateSelectionArity{expected,..}=>Some(*expected as usize),_=>None}}
#[wasm_bindgen(getter,js_name=actual)]pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionWriteError::CoordinateSelectionArity{actual,..}=>Some(*actual as usize),_=>None}}
},source=ReactionWriteErrorCause,cause=reactionwriteerror_cause);
