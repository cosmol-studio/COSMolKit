//! SMARTS and CX parser diagnostics, without string-flattened child errors.
use crate::errors::{enum_error, exception};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "StereoGroupError")]
    #[derive(Clone)]
    pub type CxStereoGroupMergeCause;
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | CxStereoGroupMergeError | PropertyUIntReadError | BondValueError | AtomPropertyError | PropertyStringError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type CxQueryLoweringCause;
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | CxQueryLoweringError | QueryGraphError | ParserCarrierError | AtomPropertyError | BondValueError | MoleculePropertyError | TemplateAttachmentOrderError"
    )]
    #[derive(Clone)]
    pub type SmartsParseCause;
}
fn merge_source(e: &ck::CxStereoGroupMergeError) -> Option<CxStereoGroupMergeCause> {
    match e {
        ck::CxStereoGroupMergeError::StereoGroup(v) => {
            Some(crate::model_errors::StereoGroupError::from_rust(*v).unchecked_into())
        }
        _ => None,
    }
}
enum_error!(CxStereoGroupMergeError,CxStereoGroupMergeErrorKind,ck::CxStereoGroupMergeError,{
    StereoGroup=>ck::CxStereoGroupMergeError::StereoGroup(_),TrackerPropertyType=>ck::CxStereoGroupMergeError::TrackerPropertyType{..},TrackerCount=>ck::CxStereoGroupMergeError::TrackerCount{..},
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<CxStereoGroupMergeCause>{merge_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn actual(&self)->Option<crate::model_errors::PropertyValueKind>{match self.inner.as_ref(){ck::CxStereoGroupMergeError::TrackerPropertyType{actual}=>Some((*actual).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn groups(&self)->Option<usize>{match self.inner.as_ref(){ck::CxStereoGroupMergeError::TrackerCount{groups,..}=>Some(*groups),_=>None}}
    #[wasm_bindgen(getter)] pub fn hashes(&self)->Option<usize>{match self.inner.as_ref(){ck::CxStereoGroupMergeError::TrackerCount{hashes,..}=>Some(*hashes),_=>None}}
},source=CxStereoGroupMergeCause,cause=merge_source);
fn cx_source(e: &ck::CxQueryLoweringError) -> Option<CxQueryLoweringCause> {
    use ck::CxQueryLoweringError as E;
    Some(
        match e {
            E::StereoGroup(v) => crate::model_errors::StereoGroupError::from_rust(*v),
            E::StereoGroupMerge(v) => CxStereoGroupMergeError::from_rust(v.clone()),
            E::Numeric(v) => crate::property_errors::PropertyUIntReadError::from_rust(v.clone()),
            E::BondProperty(v) => crate::model_errors::BondValueError::from_rust(v.clone()),
            E::AtomProperty(v) => crate::model_errors::AtomPropertyError::from_rust(v.clone()),
            E::Property(v) => crate::property_errors::PropertyStringError::from_rust(*v),
            E::MoleculeProperty(v) => {
                crate::model_errors::MoleculePropertyError::from_rust(v.clone())
            }
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(CxQueryLoweringError,CxQueryLoweringErrorKind,ck::CxQueryLoweringError,{
    StereoGroup=>ck::CxQueryLoweringError::StereoGroup(_),StereoGroupMerge=>ck::CxQueryLoweringError::StereoGroupMerge(_),Numeric=>ck::CxQueryLoweringError::Numeric(_),BondProperty=>ck::CxQueryLoweringError::BondProperty(_),AtomProperty=>ck::CxQueryLoweringError::AtomProperty(_),Property=>ck::CxQueryLoweringError::Property(_),MoleculeProperty=>ck::CxQueryLoweringError::MoleculeProperty(_),AtomIndex=>ck::CxQueryLoweringError::AtomIndex{..},BondIndex=>ck::CxQueryLoweringError::BondIndex{..},CoordinateCount=>ck::CxQueryLoweringError::CoordinateCount{..},BondAtomMismatch=>ck::CxQueryLoweringError::BondAtomMismatch{..},WedgeAtomMismatch=>ck::CxQueryLoweringError::WedgeAtomMismatch{..},UnsupportedRecord=>ck::CxQueryLoweringError::UnsupportedRecord{..},InvalidGraph=>ck::CxQueryLoweringError::InvalidGraph(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<CxQueryLoweringCause>{cx_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn index(&self)->Option<usize>{match self.inner.as_ref(){ck::CxQueryLoweringError::AtomIndex{index}|ck::CxQueryLoweringError::BondIndex{index}=>Some(*index),_=>None}}
    #[wasm_bindgen(getter)] pub fn actual(&self)->Option<usize>{match self.inner.as_ref(){ck::CxQueryLoweringError::CoordinateCount{actual,..}=>Some(*actual),_=>None}}
    #[wasm_bindgen(getter)] pub fn expected(&self)->Option<usize>{match self.inner.as_ref(){ck::CxQueryLoweringError::CoordinateCount{expected,..}=>Some(*expected),_=>None}}
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::CxQueryLoweringError::BondAtomMismatch{atom,..}|ck::CxQueryLoweringError::WedgeAtomMismatch{atom,..}=>Some(*atom),_=>None}}
    #[wasm_bindgen(getter)] pub fn bond(&self)->Option<usize>{match self.inner.as_ref(){ck::CxQueryLoweringError::BondAtomMismatch{bond,..}|ck::CxQueryLoweringError::WedgeAtomMismatch{bond,..}=>Some(*bond),_=>None}}
    #[wasm_bindgen(getter)] pub fn record(&self)->Option<String>{match self.inner.as_ref(){ck::CxQueryLoweringError::UnsupportedRecord{record}=>Some((*record).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=invalidGraph)] pub fn invalid_graph(&self)->Option<String>{match self.inner.as_ref(){ck::CxQueryLoweringError::InvalidGraph(v)=>Some(v.clone()),_=>None}}
},source=CxQueryLoweringCause,cause=cx_source);
pub(crate) fn parse_source(e: &ck::SmartsParseError) -> Option<SmartsParseCause> {
    use ck::SmartsParseError as E;
    Some(
        match e {
            E::StereoGroup(v) => crate::model_errors::StereoGroupError::from_rust(*v),
            E::CxLowering(v) => CxQueryLoweringError::from_rust(v.clone()),
            E::QueryGraph(v) => crate::query_model_errors::QueryGraphError::from_rust(v.clone()),
            E::ParserCarrier(v) => crate::parser_errors::ParserCarrierError::from_rust(v.clone()),
            E::AtomProperty(v) => crate::model_errors::AtomPropertyError::from_rust(v.clone()),
            E::BondProperty(v) => crate::model_errors::BondValueError::from_rust(v.clone()),
            E::MoleculeProperty(v) => {
                crate::model_errors::MoleculePropertyError::from_rust(v.clone())
            }
            E::TemplateAttachmentRemap { source, .. } => {
                crate::model_errors::TemplateAttachmentOrderError::from_rust(source.clone())
            }
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(SmartsParseError,SmartsParseErrorKind,ck::SmartsParseError,{
    StereoGroup=>ck::SmartsParseError::StereoGroup(_),MissingRecursiveQueryGraph=>ck::SmartsParseError::MissingRecursiveQueryGraph,CxLowering=>ck::SmartsParseError::CxLowering(_),QueryGraph=>ck::SmartsParseError::QueryGraph(_),ParserCarrier=>ck::SmartsParseError::ParserCarrier(_),AtomProperty=>ck::SmartsParseError::AtomProperty(_),BondProperty=>ck::SmartsParseError::BondProperty(_),MoleculeProperty=>ck::SmartsParseError::MoleculeProperty(_),UnclosedBracket=>ck::SmartsParseError::UnclosedBracket(_),UnexpectedCharacter=>ck::SmartsParseError::UnexpectedCharacter{..},UnexpectedEnd=>ck::SmartsParseError::UnexpectedEnd(_),InvalidAtomPrimitive=>ck::SmartsParseError::InvalidAtomPrimitive{..},UnclosedParenthesis=>ck::SmartsParseError::UnclosedParenthesis(_),UnbalancedRingClosure=>ck::SmartsParseError::UnbalancedRingClosure(_),SelfRingClosure=>ck::SmartsParseError::SelfRingClosure{..},DuplicateRingBond=>ck::SmartsParseError::DuplicateRingBond{..},CxSmiles=>ck::SmartsParseError::CxSmiles(_),Parse=>ck::SmartsParseError::Parse(_),UnsupportedFeature=>ck::SmartsParseError::UnsupportedFeature(_),TemplateAttachmentRemap=>ck::SmartsParseError::TemplateAttachmentRemap{..},
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<SmartsParseCause>{parse_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn position(&self)->Option<usize>{use ck::SmartsParseError as E;match self.inner.as_ref(){E::UnclosedBracket(v)|E::UnclosedParenthesis(v)=>Some(*v),E::UnexpectedCharacter{position,..}|E::InvalidAtomPrimitive{position,..}=>Some(*position),_=>None}}
    #[wasm_bindgen(getter)] pub fn character(&self)->Option<u8>{match self.inner.as_ref(){ck::SmartsParseError::UnexpectedCharacter{character,..}=>Some(*character),_=>None}}
    #[wasm_bindgen(getter)] pub fn context(&self)->Option<String>{match self.inner.as_ref(){ck::SmartsParseError::UnexpectedCharacter{context,..}=>Some(context.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn detail(&self)->Option<String>{match self.inner.as_ref(){ck::SmartsParseError::InvalidAtomPrimitive{detail,..}=>Some(detail.clone()),_=>None}}
    #[wasm_bindgen(getter,js_name=unexpectedEnd)] pub fn unexpected_end(&self)->Option<String>{match self.inner.as_ref(){ck::SmartsParseError::UnexpectedEnd(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter,js_name=cxSmiles)] pub fn cx_smiles(&self)->Option<String>{match self.inner.as_ref(){ck::SmartsParseError::CxSmiles(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn parse(&self)->Option<String>{match self.inner.as_ref(){ck::SmartsParseError::Parse(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter,js_name=unsupportedFeature)] pub fn unsupported_feature(&self)->Option<String>{match self.inner.as_ref(){ck::SmartsParseError::UnsupportedFeature(v)=>Some((*v).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn ring(&self)->Option<u32>{use ck::SmartsParseError as E;match self.inner.as_ref(){E::UnbalancedRingClosure(v)=>Some(*v),E::SelfRingClosure{ring,..}|E::DuplicateRingBond{ring,..}=>Some(*ring),_=>None}}
    #[wasm_bindgen(getter)] pub fn atom(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsParseError::SelfRingClosure{atom,..}=>Some(*atom),_=>None}}
    #[wasm_bindgen(getter,js_name=beginAtom)] pub fn begin_atom(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsParseError::DuplicateRingBond{begin_atom,..}=>Some(*begin_atom),_=>None}}
    #[wasm_bindgen(getter,js_name=endAtom)] pub fn end_atom(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsParseError::DuplicateRingBond{end_atom,..}=>Some(*end_atom),_=>None}}
    #[wasm_bindgen(getter)] pub fn carrier(&self)->Option<usize>{match self.inner.as_ref(){ck::SmartsParseError::TemplateAttachmentRemap{carrier,..}=>Some(*carrier),_=>None}}
},source=SmartsParseCause,cause=parse_source);
