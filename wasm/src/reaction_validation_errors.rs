use crate::model_errors::AtomPropertyError;
use crate::property_errors::PropertyIntReadError;
use crate::{
    errors::{enum_error, exception},
    reaction_validation_values::{ReactionRole, ReactionValidationReport},
    storage::ReadonlyValue,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
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
