//! Pattern fingerprint facade and its exact domain error.
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fingerprint_values::Fingerprint,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen(js_name=fingerprintPatternQuery)]
pub fn fingerprint_pattern_query(
    query: crate::query_values::QueryGraphInput,
    params: Option<PatternFingerprintParamsInput>,
) -> Result<Fingerprint, JsValue> {
    crate::query_values::QueryGraph::with_ref(&query, |query| {
        PatternFingerprintParams::with_input(params.as_ref(), |p| {
            ck::fingerprint_pattern_query_with_params(query, p)
                .map(Fingerprint)
                .map_err(PatternFingerprintError::from_rust)
        })
    })
}
#[wasm_bindgen(js_name=fingerprintPatternQueryWithParams)]
pub fn fingerprint_pattern_query_with_params(
    query: crate::query_values::QueryGraphInput,
    params: PatternFingerprintParamsInput,
) -> Result<Fingerprint, JsValue> {
    fingerprint_pattern_query(query, Some(params))
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | QueryGraphError | QueryAtomConversionError | RingFindingError | SmartsParseError | QueryCompileError | QueryMatchContextError | SubstructMatchError | FingerprintError"
    )]
    #[derive(Clone)]
    pub type PatternCause;
}
fn source(e: &ck::PatternFingerprintError) -> Option<PatternCause> {
    use ck::PatternFingerprintError as E;
    Some(
        match e {
            E::Topology(v) => crate::model_errors::TopologyValidationError::from_rust(v.clone()),
            E::Query(v) => crate::query_model_errors::QueryGraphError::from_rust(v.clone()),
            E::QueryCarrier(v) => {
                crate::query_model_errors::QueryAtomConversionError::from_rust(v.clone())
            }
            E::Rings(v) => crate::ring_errors::RingFindingError::from_rust(v.clone()),
            E::Smarts(v) => crate::smarts_parse_errors::SmartsParseError::from_rust(v.clone()),
            E::QueryCompile(v) => {
                crate::query_compile_error::QueryCompileError::from_rust(v.clone())
            }
            E::QueryContext(v) => crate::match_errors::QueryMatchContextError::from_rust(v.clone()),
            E::Match(v) => crate::match_errors::SubstructMatchError::from_rust(v.clone()),
            E::Value(v) => crate::fingerprint_errors::FingerprintError::from_rust(v.clone()),
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(PatternFingerprintError,PatternFingerprintErrorKind,ck::PatternFingerprintError,{
    EmptyFingerprint=>ck::PatternFingerprintError::EmptyFingerprint,InvalidArguments=>ck::PatternFingerprintError::InvalidArguments{..},BitLengthMismatch=>ck::PatternFingerprintError::BitLengthMismatch{..},Topology=>ck::PatternFingerprintError::Topology(_),Query=>ck::PatternFingerprintError::Query(_),QueryCarrier=>ck::PatternFingerprintError::QueryCarrier(_),Rings=>ck::PatternFingerprintError::Rings(_),Smarts=>ck::PatternFingerprintError::Smarts(_),QueryCompile=>ck::PatternFingerprintError::QueryCompile(_),QueryContext=>ck::PatternFingerprintError::QueryContext(_),Match=>ck::PatternFingerprintError::Match(_),Value=>ck::PatternFingerprintError::Value(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<PatternCause>{source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::PatternFingerprintError::InvalidArguments{reason}=>Some((*reason).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn left(&self)->Option<usize>{match self.inner.as_ref(){ck::PatternFingerprintError::BitLengthMismatch{left,..}=>Some(*left),_=>None}}
    #[wasm_bindgen(getter)] pub fn right(&self)->Option<usize>{match self.inner.as_ref(){ck::PatternFingerprintError::BitLengthMismatch{right,..}=>Some(*right),_=>None}}
},source=PatternCause,cause=source);
ck::__ck_configuration_PatternFingerprintParams!(
    crate::project_parameters,
    ck::PatternFingerprintParams,
    {},
    { PatternFingerprintParams }
);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintPattern)]
    pub fn fingerprint_pattern(
        &self,
        params: Option<PatternFingerprintParamsInput>,
    ) -> Result<Fingerprint, JsValue> {
        PatternFingerprintParams::with_input(params.as_ref(), |params| {
            self.0
                .fingerprint_pattern_with_params(params)
                .map(Fingerprint)
                .map_err(PatternFingerprintError::from_rust)
        })
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = fingerprintPatternWithParams)]
    pub fn fingerprint_pattern_with_params(
        &self,
        params: PatternFingerprintParamsInput,
    ) -> Result<Fingerprint, JsValue> {
        self.fingerprint_pattern(Some(params))
    }
}
