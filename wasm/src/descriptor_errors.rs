//! Descriptor errors preserve prepared-state failures and actual algorithm sources.
use crate::errors::{enum_error, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "SmartsParseError | SubstructMatchError | QueryMatchContextError"
    )]
    #[derive(Clone)]
    pub type DescriptorSearchSource;
    #[wasm_bindgen(
        typescript_type = "TopologyEditError | TopologyValidationError | PathError | HydrogenError | ValenceError | RingFindingError | DescriptorSearchCause | LegacyStereoError"
    )]
    #[derive(Clone)]
    pub type DescriptorSource;
    #[wasm_bindgen(typescript_type = "DescriptorError")]
    #[derive(Clone)]
    pub type DescriptorReadSource;
}
fn search_source(
    owner: &ReadonlyValue<ck::DescriptorSearchCause>,
) -> Option<DescriptorSearchSource> {
    use crate::errors::child;
    use ck::DescriptorSearchCause as E;
    Some(match owner.as_ref() {
        E::Compile(_) => child!(
            owner,
            E::Compile,
            crate::smarts_parse_errors::SmartsParseError
        ),
        E::Match(_) => child!(owner, E::Match, crate::match_errors::SubstructMatchError),
        E::Context(_) => child!(
            owner,
            E::Context,
            crate::match_errors::QueryMatchContextError
        ),
    })
}
enum_error!(DescriptorSearchCause,DescriptorSearchCauseKind,ck::DescriptorSearchCause,{
    Compile=>ck::DescriptorSearchCause::Compile(_),Match=>ck::DescriptorSearchCause::Match(_),Context=>ck::DescriptorSearchCause::Context(_)
},{#[wasm_bindgen(getter)]pub fn source(&self)->Option<DescriptorSearchSource>{search_source(&self.inner)}},source=DescriptorSearchSource,cause=search_source);
macro_rules! source_field {
    ($owner:expr,$variant:path,$binding:path) => {{
        exception(<$binding>::from_owner($owner.child(|e| match e {
            $variant { source, .. } => source,
            _ => unreachable!("checked immutable native variant"),
        })))
        .unchecked_into()
    }};
}
fn source(owner: &ReadonlyValue<ck::DescriptorError>) -> Option<DescriptorSource> {
    use ck::DescriptorError as E;
    Some(match owner.as_ref() {
        E::TopologyEdit { .. } => source_field!(
            owner,
            E::TopologyEdit,
            crate::topology_errors::TopologyEditError
        ),
        E::InvalidTopology { .. } => source_field!(
            owner,
            E::InvalidTopology,
            crate::model_errors::TopologyValidationError
        ),
        E::Path { .. } => source_field!(owner, E::Path, crate::path_errors::PathError),
        E::Hydrogens { .. } => {
            source_field!(owner, E::Hydrogens, crate::hydrogen_errors::HydrogenError)
        }
        E::Valence { .. } => source_field!(owner, E::Valence, crate::valence_errors::ValenceError),
        E::Ring { .. } => source_field!(owner, E::Ring, crate::ring_errors::RingFindingError),
        E::Search { .. } => source_field!(owner, E::Search, DescriptorSearchCause),
        E::Stereo { .. } => source_field!(
            owner,
            E::Stereo,
            crate::legacy_stereo_errors::LegacyStereoError
        ),
        _ => return None,
    })
}
enum_error!(DescriptorError,DescriptorErrorKind,ck::DescriptorError,{
    TopologyEdit=>ck::DescriptorError::TopologyEdit{..},MissingFinalHydrogenState=>ck::DescriptorError::MissingFinalHydrogenState{..},InvalidConnectivityPath=>ck::DescriptorError::InvalidConnectivityPath{..},InvalidTopology=>ck::DescriptorError::InvalidTopology{..},Path=>ck::DescriptorError::Path{..},InvalidHallKierContributionRows=>ck::DescriptorError::InvalidHallKierContributionRows{..},MissingComputedScalar=>ck::DescriptorError::MissingComputedScalar{..},MissingLabuteHydrogens=>ck::DescriptorError::MissingLabuteHydrogens{..},MissingLabuteAsa=>ck::DescriptorError::MissingLabuteAsa{..},MismatchedBinArrays=>ck::DescriptorError::MismatchedBinArrays{..},CrippenParamNumeric=>ck::DescriptorError::CrippenParamNumeric{..},MissingCrippenMrContributions=>ck::DescriptorError::MissingCrippenMrContributions{..},MissingCrippenMr=>ck::DescriptorError::MissingCrippenMr{..},InvalidCrippenOptionalRows=>ck::DescriptorError::InvalidCrippenOptionalRows{..},MissingCrippenDefaultPattern=>ck::DescriptorError::MissingCrippenDefaultPattern{..},Hydrogens=>ck::DescriptorError::Hydrogens{..},CountOverflow=>ck::DescriptorError::CountOverflow{..},Valence=>ck::DescriptorError::Valence{..},InvalidValenceRows=>ck::DescriptorError::InvalidValenceRows{..},Unsupported=>ck::DescriptorError::Unsupported{..},Ring=>ck::DescriptorError::Ring{..},Search=>ck::DescriptorError::Search{..},Stereo=>ck::DescriptorError::Stereo{..}
},{
    #[wasm_bindgen(getter)]pub fn source(&self)->Option<DescriptorSource>{source(&self.inner)}
    #[wasm_bindgen(getter)]pub fn function(&self)->Option<String>{use ck::DescriptorError as E;match self.inner.as_ref(){E::TopologyEdit{function,..}|E::InvalidConnectivityPath{function,..}|E::InvalidTopology{function,..}|E::Path{function,..}|E::MissingComputedScalar{function,..}|E::MissingLabuteHydrogens{function}|E::MissingLabuteAsa{function}|E::MissingCrippenMrContributions{function}|E::MissingCrippenMr{function}|E::Hydrogens{function,..}|E::CountOverflow{function,..}|E::Valence{function,..}|E::InvalidValenceRows{function,..}|E::Unsupported{function,..}|E::Ring{function,..}|E::Search{function,..}|E::Stereo{function,..}=>Some((*function).into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn field(&self)->Option<String>{use ck::DescriptorError as E;match self.inner.as_ref(){E::MissingFinalHydrogenState{field}|E::CrippenParamNumeric{field,..}|E::InvalidCrippenOptionalRows{field,..}|E::CountOverflow{field,..}|E::InvalidValenceRows{field,..}=>Some((*field).into()),_=>None}}
    #[wasm_bindgen(getter)]pub fn actual(&self)->Option<usize>{use ck::DescriptorError as E;match self.inner.as_ref(){E::InvalidHallKierContributionRows{actual,..}|E::InvalidCrippenOptionalRows{actual,..}|E::InvalidValenceRows{actual,..}=>Some(*actual),_=>None}}
    #[wasm_bindgen(getter)]pub fn expected(&self)->Option<usize>{use ck::DescriptorError as E;match self.inner.as_ref(){E::InvalidCrippenOptionalRows{expected,..}|E::InvalidValenceRows{expected,..}=>Some(*expected),_=>None}}
    #[wasm_bindgen(getter)]pub fn minimum(&self)->Option<usize>{match self.inner.as_ref(){ck::DescriptorError::InvalidHallKierContributionRows{minimum,..}=>Some(*minimum),_=>None}}
    #[wasm_bindgen(getter,js_name=expectedRows)]pub fn expected_rows(&self)->Option<usize>{match self.inner.as_ref(){ck::DescriptorError::InvalidConnectivityPath{expected_rows,..}=>Some(*expected_rows),_=>None}}
    #[wasm_bindgen(getter,js_name=actualRows)]pub fn actual_rows(&self)->Option<usize>{match self.inner.as_ref(){ck::DescriptorError::InvalidConnectivityPath{actual_rows,..}=>*actual_rows,_=>None}}
    #[wasm_bindgen(getter,js_name=includeSulfurPhosphorus)]pub fn include_sulfur_phosphorus(&self)->Option<bool>{match self.inner.as_ref(){ck::DescriptorError::MissingComputedScalar{include_sulfur_phosphorus,..}=>Some(*include_sulfur_phosphorus),_=>None}}
    #[wasm_bindgen(getter,js_name=contribsLen)]pub fn contribs_len(&self)->Option<usize>{match self.inner.as_ref(){ck::DescriptorError::MismatchedBinArrays{contribs_len,..}=>Some(*contribs_len),_=>None}}
    #[wasm_bindgen(getter,js_name=binPropLen)]pub fn bin_prop_len(&self)->Option<usize>{match self.inner.as_ref(){ck::DescriptorError::MismatchedBinArrays{bin_prop_len,..}=>Some(*bin_prop_len),_=>None}}
    #[wasm_bindgen(getter,js_name=binsLen)]pub fn bins_len(&self)->Option<usize>{match self.inner.as_ref(){ck::DescriptorError::MismatchedBinArrays{bins_len,..}=>Some(*bins_len),_=>None}}
    #[wasm_bindgen(getter)]pub fn cell(&self)->Option<String>{match self.inner.as_ref(){ck::DescriptorError::CrippenParamNumeric{cell,..}=>Some(cell.clone()),_=>None}}
    #[wasm_bindgen(getter)]pub fn row(&self)->Option<u32>{match self.inner.as_ref(){ck::DescriptorError::MissingCrippenDefaultPattern{row}=>Some(*row),_=>None}}
    #[wasm_bindgen(getter)]pub fn detail(&self)->Option<String>{match self.inner.as_ref(){ck::DescriptorError::Unsupported{detail,..}=>Some(detail.clone()),_=>None}}
},source=DescriptorSource,cause=source);
fn read_source(owner: &ReadonlyValue<ck::DescriptorReadError>) -> Option<DescriptorReadSource> {
    match owner.as_ref() {
        ck::DescriptorReadError::Algorithm { .. } => Some(source_field!(
            owner,
            ck::DescriptorReadError::Algorithm,
            DescriptorError
        )),
        _ => None,
    }
}
enum_error!(DescriptorReadError,DescriptorReadErrorKind,ck::DescriptorReadError,{
    MissingPreparedValence=>ck::DescriptorReadError::MissingPreparedValence,CachePoisoned=>ck::DescriptorReadError::CachePoisoned,MissingInitializedRings=>ck::DescriptorReadError::MissingInitializedRings,Algorithm=>ck::DescriptorReadError::Algorithm{..}
},{#[wasm_bindgen(getter)]pub fn source(&self)->Option<DescriptorReadSource>{read_source(&self.inner)}},source=DescriptorReadSource,cause=read_source);
