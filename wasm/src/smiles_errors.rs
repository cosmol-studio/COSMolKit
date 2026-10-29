//! Public SMILES construction and serialization errors preserve native layers.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | ValenceError | RingFindingError | LegacyStereoError | CanonicalRankError | SmilesParseError | PropertyStringError | AtomPropertyError | BondValueError"
    )]
    #[derive(Clone)]
    pub type FragmentWriteInputErrorCause;
    #[wasm_bindgen(typescript_type = "StereoGroupError")]
    pub type FragmentWriteInputErrorStereoGroupPayload;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type FragmentWriteInputErrorValencePayload;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type FragmentWriteInputErrorRingPayload;
    #[wasm_bindgen(typescript_type = "LegacyStereoError")]
    pub type FragmentWriteInputErrorStereoPayload;
    #[wasm_bindgen(typescript_type = "CanonicalRankError")]
    pub type FragmentWriteInputErrorRankPayload;
    #[wasm_bindgen(typescript_type = "SmilesParseError")]
    pub type FragmentWriteInputErrorWriterPayload;
    #[wasm_bindgen(typescript_type = "PropertyStringError")]
    pub type FragmentWriteInputErrorCipPropertyPayload;
    #[wasm_bindgen(typescript_type = "AtomPropertyError")]
    pub type FragmentWriteInputErrorAtomPropertyPayload;
    #[wasm_bindgen(typescript_type = "BondValueError")]
    pub type FragmentWriteInputErrorStereoBondPayload;
}
fn fragment_write_input_error_cause(
    error: &ck::FragmentWriteInputError,
) -> Option<FragmentWriteInputErrorCause> {
    match error {
        ck::FragmentWriteInputError::StereoGroup(v) => {
            Some(crate::model_errors::StereoGroupError::from_rust(v.clone()).unchecked_into())
        }
        ck::FragmentWriteInputError::Valence(v) => {
            Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into())
        }
        ck::FragmentWriteInputError::Ring(v) => {
            Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into())
        }
        ck::FragmentWriteInputError::Stereo(v) => Some(
            crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::FragmentWriteInputError::Rank(v) => {
            Some(crate::canonical_errors::CanonicalRankError::from_rust(v.clone()).unchecked_into())
        }
        ck::FragmentWriteInputError::Writer(v) => Some(
            crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::FragmentWriteInputError::CipProperty(v) => {
            Some(crate::property_errors::PropertyStringError::from_rust(v.clone()).unchecked_into())
        }
        ck::FragmentWriteInputError::AtomProperty(v) => {
            Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into())
        }
        ck::FragmentWriteInputError::StereoBond(v) => {
            Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(FragmentWriteInputError,FragmentWriteInputErrorKind,ck::FragmentWriteInputError,{
StereoGroup=>ck::FragmentWriteInputError::StereoGroup(_),
Valence=>ck::FragmentWriteInputError::Valence(_),
Ring=>ck::FragmentWriteInputError::Ring(_),
Stereo=>ck::FragmentWriteInputError::Stereo(_),
Rank=>ck::FragmentWriteInputError::Rank(_),
Writer=>ck::FragmentWriteInputError::Writer(_),
CipProperty=>ck::FragmentWriteInputError::CipProperty(_),
AtomProperty=>ck::FragmentWriteInputError::AtomProperty(_),
StereoBond=>ck::FragmentWriteInputError::StereoBond(_),
RootAtomOutOfRange=>ck::FragmentWriteInputError::RootAtomOutOfRange{..},
RootAtomNotSelected=>ck::FragmentWriteInputError::RootAtomNotSelected{..},
AtomSymbolsTooShort=>ck::FragmentWriteInputError::AtomSymbolsTooShort{..},
BondSymbolsTooShort=>ck::FragmentWriteInputError::BondSymbolsTooShort{..},
AtomOutOfRange=>ck::FragmentWriteInputError::AtomOutOfRange{..},
BondOutOfRange=>ck::FragmentWriteInputError::BondOutOfRange{..},
RootAtomRequiresSingleFragment=>ck::FragmentWriteInputError::RootAtomRequiresSingleFragment{..},
ValenceAssignmentLengthMismatch=>ck::FragmentWriteInputError::ValenceAssignmentLengthMismatch{..},
NoAtomsProvided=>ck::FragmentWriteInputError::NoAtomsProvided,
InvalidTopology=>ck::FragmentWriteInputError::InvalidTopology(_),
},{
#[wasm_bindgen(getter,js_name=stereoGroup)] pub fn stereo_group(&self)->Option<FragmentWriteInputErrorStereoGroupPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::StereoGroup(v)=>Some(crate::model_errors::StereoGroupError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=valence)] pub fn valence(&self)->Option<FragmentWriteInputErrorValencePayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::Valence(v)=>Some(crate::valence_errors::ValenceError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=ring)] pub fn ring(&self)->Option<FragmentWriteInputErrorRingPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::Ring(v)=>Some(crate::ring_errors::RingFindingError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stereo)] pub fn stereo(&self)->Option<FragmentWriteInputErrorStereoPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::Stereo(v)=>Some(crate::legacy_stereo_errors::LegacyStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=rank)] pub fn rank(&self)->Option<FragmentWriteInputErrorRankPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::Rank(v)=>Some(crate::canonical_errors::CanonicalRankError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=writer)] pub fn writer(&self)->Option<FragmentWriteInputErrorWriterPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::Writer(v)=>Some(crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=cipProperty)] pub fn cip_property(&self)->Option<FragmentWriteInputErrorCipPropertyPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::CipProperty(v)=>Some(crate::property_errors::PropertyStringError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atomProperty)] pub fn atom_property(&self)->Option<FragmentWriteInputErrorAtomPropertyPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::AtomProperty(v)=>Some(crate::model_errors::AtomPropertyError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stereoBond)] pub fn stereo_bond(&self)->Option<FragmentWriteInputErrorStereoBondPayload> {match self.inner.as_ref(){ck::FragmentWriteInputError::StereoBond(v)=>Some(crate::model_errors::BondValueError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=atomIndex)] pub fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::RootAtomOutOfRange{atom_index,..}=>Some(*atom_index),
ck::FragmentWriteInputError::RootAtomNotSelected{atom_index,..}=>Some(*atom_index),
ck::FragmentWriteInputError::AtomOutOfRange{atom_index,..}=>Some(*atom_index),
ck::FragmentWriteInputError::RootAtomRequiresSingleFragment{atom_index,..}=>Some(*atom_index),
_=>None}}
#[wasm_bindgen(getter,js_name=atomCount)] pub fn atom_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::RootAtomOutOfRange{atom_count,..}=>Some(*atom_count),
ck::FragmentWriteInputError::AtomOutOfRange{atom_count,..}=>Some(*atom_count),
ck::FragmentWriteInputError::ValenceAssignmentLengthMismatch{atom_count,..}=>Some(*atom_count),
_=>None}}
#[wasm_bindgen(getter,js_name=len)] pub fn len(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::AtomSymbolsTooShort{len,..}=>Some(*len),
ck::FragmentWriteInputError::BondSymbolsTooShort{len,..}=>Some(*len),
_=>None}}
#[wasm_bindgen(getter,js_name=expected)] pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::AtomSymbolsTooShort{expected,..}=>Some(*expected),
ck::FragmentWriteInputError::BondSymbolsTooShort{expected,..}=>Some(*expected),
_=>None}}
#[wasm_bindgen(getter,js_name=bondIndex)] pub fn bond_index(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::BondOutOfRange{bond_index,..}=>Some(*bond_index),
_=>None}}
#[wasm_bindgen(getter,js_name=bondCount)] pub fn bond_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::BondOutOfRange{bond_count,..}=>Some(*bond_count),
_=>None}}
#[wasm_bindgen(getter,js_name=fragmentCount)] pub fn fragment_count(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::RootAtomRequiresSingleFragment{fragment_count,..}=>Some(*fragment_count),
_=>None}}
#[wasm_bindgen(getter,js_name=explicitLen)] pub fn explicit_len(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::ValenceAssignmentLengthMismatch{explicit_len,..}=>Some(*explicit_len),
_=>None}}
#[wasm_bindgen(getter,js_name=implicitLen)] pub fn implicit_len(&self)->Option<usize>{match self.inner.as_ref(){
ck::FragmentWriteInputError::ValenceAssignmentLengthMismatch{implicit_len,..}=>Some(*implicit_len),
_=>None}}
#[wasm_bindgen(getter,js_name=invalidTopology)] pub fn invalid_topology(&self)->Option<String> {match self.inner.as_ref(){ck::FragmentWriteInputError::InvalidTopology(v)=>Some(v.clone()),_=>None}}
}, source=FragmentWriteInputErrorCause, cause=fragment_write_input_error_cause);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "SmilesParseError | HydrogenError | SanitizeError | SmilesStereoError | OperationError"
    )]
    #[derive(Clone)]
    pub type SmilesErrorCause;
    #[wasm_bindgen(typescript_type = "SmilesParseError")]
    pub type SmilesErrorParsePayload;
    #[wasm_bindgen(typescript_type = "HydrogenError")]
    pub type SmilesErrorHydrogenPayload;
    #[wasm_bindgen(typescript_type = "SanitizeError")]
    pub type SmilesErrorSanitizePayload;
    #[wasm_bindgen(typescript_type = "SmilesStereoError")]
    pub type SmilesErrorStereoPayload;
    #[wasm_bindgen(typescript_type = "OperationError")]
    pub type SmilesErrorConstructionPayload;
}
fn smiles_error_cause(error: &ck::SmilesError) -> Option<SmilesErrorCause> {
    match error {
        ck::SmilesError::Parse(v) => Some(
            crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesError::Hydrogen(v) => {
            Some(crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesError::Sanitize(v) => {
            Some(crate::sanitize_errors::SanitizeError::from_rust(v.clone()).unchecked_into())
        }
        ck::SmilesError::Stereo(v) => Some(
            crate::smiles_stereo_errors::SmilesStereoError::from_rust(v.clone()).unchecked_into(),
        ),
        ck::SmilesError::Construction(v) => {
            Some(crate::operation_errors::OperationError::from_rust(v.clone()).unchecked_into())
        }
    }
}
enum_error!(SmilesError,SmilesErrorKind,ck::SmilesError,{
Parse=>ck::SmilesError::Parse(_),
Hydrogen=>ck::SmilesError::Hydrogen(_),
Sanitize=>ck::SmilesError::Sanitize(_),
Stereo=>ck::SmilesError::Stereo(_),
Construction=>ck::SmilesError::Construction(_),
},{
#[wasm_bindgen(getter,js_name=parse)] pub fn parse(&self)->Option<SmilesErrorParsePayload> {match self.inner.as_ref(){ck::SmilesError::Parse(v)=>Some(crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=hydrogen)] pub fn hydrogen(&self)->Option<SmilesErrorHydrogenPayload> {match self.inner.as_ref(){ck::SmilesError::Hydrogen(v)=>Some(crate::hydrogen_errors::HydrogenError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=sanitize)] pub fn sanitize(&self)->Option<SmilesErrorSanitizePayload> {match self.inner.as_ref(){ck::SmilesError::Sanitize(v)=>Some(crate::sanitize_errors::SanitizeError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=stereo)] pub fn stereo(&self)->Option<SmilesErrorStereoPayload> {match self.inner.as_ref(){ck::SmilesError::Stereo(v)=>Some(crate::smiles_stereo_errors::SmilesStereoError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=construction)] pub fn construction(&self)->Option<SmilesErrorConstructionPayload> {match self.inner.as_ref(){ck::SmilesError::Construction(v)=>Some(crate::operation_errors::OperationError::from_rust(v.clone()).unchecked_into()),_=>None}}
}, source=SmilesErrorCause, cause=smiles_error_cause);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "SmilesParseError | FragmentWriteInputError")]
    #[derive(Clone)]
    pub type SmilesWriteErrorCause;
    #[wasm_bindgen(typescript_type = "SmilesParseError")]
    pub type SmilesWriteErrorWritePayload;
    #[wasm_bindgen(typescript_type = "FragmentWriteInputError")]
    pub type SmilesWriteErrorFragmentPayload;
}
fn smiles_write_error_cause(
    error: &crate::storage::ReadonlyValue<ck::SmilesWriteError>,
) -> Option<SmilesWriteErrorCause> {
    match error.as_ref() {
        ck::SmilesWriteError::Write(_) => Some(
            crate::smiles_parse_errors::SmilesParseError::from_owner(error.child(
                |root| match root {
                    ck::SmilesWriteError::Write(value) => value,
                    _ => unreachable!("immutable Write variant selected"),
                },
            ))
            .unchecked_into(),
        ),
        ck::SmilesWriteError::Fragment(_) => Some(
            crate::smiles_errors::FragmentWriteInputError::from_owner(error.child(
                |root| match root {
                    ck::SmilesWriteError::Fragment(value) => value,
                    _ => unreachable!("immutable Fragment variant selected"),
                },
            ))
            .unchecked_into(),
        ),
    }
}
enum_error!(SmilesWriteError,SmilesWriteErrorKind,ck::SmilesWriteError,{
Write=>ck::SmilesWriteError::Write(_),
Fragment=>ck::SmilesWriteError::Fragment(_),
},{
#[wasm_bindgen(getter,js_name=write)] pub fn write(&self)->Option<SmilesWriteErrorWritePayload> {match self.inner.as_ref(){ck::SmilesWriteError::Write(v)=>Some(crate::smiles_parse_errors::SmilesParseError::from_rust(v.clone()).unchecked_into()),_=>None}}
#[wasm_bindgen(getter,js_name=fragment)] pub fn fragment(&self)->Option<SmilesWriteErrorFragmentPayload> {match self.inner.as_ref(){ck::SmilesWriteError::Fragment(v)=>Some(crate::smiles_errors::FragmentWriteInputError::from_rust(v.clone()).unchecked_into()),_=>None}}
}, source=SmilesWriteErrorCause, cause=smiles_write_error_cause);
