//! MACCS calls share canonical configuration and exact native diagnostics.
use crate::query_compile_error::QueryCompileError;
use crate::{
    Molecule,
    errors::{enum_error, exception},
    fingerprint_values::Fingerprint,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "TopologyValidationError | RingFindingError | ValenceError | PathError | SmartsParseError | QueryCompileError | QueryMatchContextError | SubstructMatchError | FingerprintError"
    )]
    #[derive(Clone)]
    pub type MaccsCause;
}
fn source(e: &ck::MaccsFingerprintError) -> Option<MaccsCause> {
    use ck::MaccsFingerprintError as E;
    Some(
        match e {
            E::Topology(v) => crate::model_errors::TopologyValidationError::from_rust(v.clone()),
            E::Rings(v) => crate::ring_errors::RingFindingError::from_rust(v.clone()),
            E::Valence(v) => crate::valence_errors::ValenceError::from_rust(v.clone()),
            E::Paths(v) => crate::path_errors::PathError::from_rust(v.clone()),
            E::Smarts(v) => crate::smarts_parse_errors::SmartsParseError::from_rust(v.clone()),
            E::QueryCompile(v) => QueryCompileError::from_rust(v.clone()),
            E::QueryContext(v) => crate::match_errors::QueryMatchContextError::from_rust(v.clone()),
            E::Match(v) => crate::match_errors::SubstructMatchError::from_rust(v.clone()),
            E::Value(v) => crate::fingerprint_errors::FingerprintError::from_rust(v.clone()),
            _ => return None,
        }
        .unchecked_into(),
    )
}
enum_error!(MaccsFingerprintError,MaccsFingerprintErrorKind,ck::MaccsFingerprintError,{
    UnsupportedOption=>ck::MaccsFingerprintError::UnsupportedOption{..},MissingPattern=>ck::MaccsFingerprintError::MissingPattern{..},Topology=>ck::MaccsFingerprintError::Topology(_),Rings=>ck::MaccsFingerprintError::Rings(_),Valence=>ck::MaccsFingerprintError::Valence(_),Paths=>ck::MaccsFingerprintError::Paths(_),Smarts=>ck::MaccsFingerprintError::Smarts(_),QueryCompile=>ck::MaccsFingerprintError::QueryCompile(_),QueryContext=>ck::MaccsFingerprintError::QueryContext(_),Match=>ck::MaccsFingerprintError::Match(_),Value=>ck::MaccsFingerprintError::Value(_),
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<MaccsCause>{source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn option(&self)->Option<String>{match self.inner.as_ref(){ck::MaccsFingerprintError::UnsupportedOption{option,..}=>Some((*option).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::MaccsFingerprintError::UnsupportedOption{reason,..}=>Some((*reason).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn bit(&self)->Option<usize>{match self.inner.as_ref(){ck::MaccsFingerprintError::MissingPattern{bit}=>Some(*bit),_=>None}}
},source=MaccsCause,cause=source);
ck::__ck_configuration_MaccsFingerprintParams!(
    crate::project_parameters,
    ck::MaccsFingerprintParams,
    {},
    { MaccsFingerprintParams }
);
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintMaccs)]
    pub fn fingerprint_maccs(
        &self,
        params: Option<MaccsFingerprintParamsInput>,
    ) -> Result<Fingerprint, JsValue> {
        MaccsFingerprintParams::with_input(params.as_ref(), |params| {
            self.0
                .fingerprint_maccs_with_params(params)
                .map(Fingerprint)
                .map_err(MaccsFingerprintError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=fingerprintMaccsRaw)]
    pub fn fingerprint_maccs_raw(&self) -> Result<Fingerprint, JsValue> {
        self.0
            .fingerprint_maccs_raw()
            .map(Fingerprint)
            .map_err(MaccsFingerprintError::from_rust)
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = fingerprintMaccsWithParams)]
    pub fn fingerprint_maccs_with_params(
        &self,
        params: MaccsFingerprintParamsInput,
    ) -> Result<Fingerprint, JsValue> {
        self.fingerprint_maccs(Some(params))
    }
}
