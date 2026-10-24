//! Native fingerprint and atom-code errors preserve their numeric payload widths.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ValenceError | CipLabelerError | PropertyStringError")]
    #[derive(Clone)]
    pub type AtomCodeCause;
}
fn atom_code_source(e: &ck::AtomCodeError) -> Option<AtomCodeCause> {
    use ck::AtomCodeError as E;
    Some(
        match e {
            E::Valence(v) => crate::valence_errors::ValenceError::from_rust(v.clone()),
            E::Cip(v) => crate::cip_labeler_errors::CipLabelerError::from_rust(v.clone()),
            E::PropertyString(v) => crate::property_errors::PropertyStringError::from_rust(*v),
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(AtomCodeError,AtomCodeErrorKind,ck::AtomCodeError,{
    Precondition=>ck::AtomCodeError::Precondition{..},Postcondition=>ck::AtomCodeError::Postcondition{..},Valence=>ck::AtomCodeError::Valence(_),Cip=>ck::AtomCodeError::Cip(_),PropertyString=>ck::AtomCodeError::PropertyString(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<AtomCodeCause>{atom_code_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn what(&self)->Option<String>{match self.inner.as_ref(){ck::AtomCodeError::Precondition{what}|ck::AtomCodeError::Postcondition{what,..}=>Some((*what).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn code(&self)->Option<u32>{match self.inner.as_ref(){ck::AtomCodeError::Postcondition{code,..}=>Some(*code),_=>None}}
},source=AtomCodeCause,cause=atom_code_source);
enum_error!(FingerprintError,FingerprintErrorKind,ck::FingerprintError,{
    EmptyFingerprint=>ck::FingerprintError::EmptyFingerprint,Unsupported=>ck::FingerprintError::Unsupported,SparseIndexOutOfRange=>ck::FingerprintError::SparseIndexOutOfRange{..},BitLengthMismatch=>ck::FingerprintError::BitLengthMismatch{..},InvalidFoldFactor=>ck::FingerprintError::InvalidFoldFactor{..},RangeError=>ck::FingerprintError::RangeError{..},UndefinedArithmetic=>ck::FingerprintError::UndefinedArithmetic{..},PreconditionViolation=>ck::FingerprintError::PreconditionViolation{..},InvalidArguments=>ck::FingerprintError::InvalidArguments{..},
},{
    #[wasm_bindgen(getter)] pub fn index(&self)->Option<u64>{match self.inner.as_ref(){ck::FingerprintError::SparseIndexOutOfRange{index,..}=>Some(*index),_=>None}}
    #[wasm_bindgen(getter)] pub fn size(&self)->Option<u64>{match self.inner.as_ref(){ck::FingerprintError::SparseIndexOutOfRange{size,..}=>Some(*size),_=>None}}
    #[wasm_bindgen(getter)] pub fn left(&self)->Option<u64>{match self.inner.as_ref(){ck::FingerprintError::BitLengthMismatch{left,..}=>Some(*left),_=>None}}
    #[wasm_bindgen(getter)] pub fn right(&self)->Option<u64>{match self.inner.as_ref(){ck::FingerprintError::BitLengthMismatch{right,..}=>Some(*right),_=>None}}
    #[wasm_bindgen(getter)] pub fn factor(&self)->Option<u32>{match self.inner.as_ref(){ck::FingerprintError::InvalidFoldFactor{factor,..}=>Some(*factor),_=>None}}
    #[wasm_bindgen(getter,js_name=nBits)] pub fn n_bits(&self)->Option<u32>{match self.inner.as_ref(){ck::FingerprintError::InvalidFoldFactor{n_bits,..}=>Some(*n_bits),_=>None}}
    #[wasm_bindgen(getter)] pub fn value(&self)->Option<f64>{match self.inner.as_ref(){ck::FingerprintError::RangeError{value}=>Some(*value),_=>None}}
    #[wasm_bindgen(getter)] pub fn site(&self)->Option<String>{match self.inner.as_ref(){ck::FingerprintError::UndefinedArithmetic{site}=>Some((*site).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn what(&self)->Option<String>{match self.inner.as_ref(){ck::FingerprintError::PreconditionViolation{what}=>Some((*what).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::FingerprintError::InvalidArguments{reason}=>Some((*reason).into()),_=>None}}
});
