//! Exact native fingerprint error graph, retaining non-Clone JSON causes.
use crate::errors::{enum_error, error_projection, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
pub enum SerdeJsonCategory {
    Io,
    Syntax,
    Data,
    Eof,
}
fn json_cause(
    owner: &ReadonlyValue<serde_json::Error>,
) -> Option<crate::external_errors::StandardCause> {
    owner.source().map(|_| {
        crate::external_errors::standard_source(
            owner.child(|e| e.source().expect("existing immutable JSON source")),
        )
        .unchecked_into()
    })
}
error_projection!(
    SerdeJsonError,
    serde_json::Error,
    {
        #[wasm_bindgen(getter)]
        pub fn line(&self) -> usize {
            self.inner.line()
        }
        #[wasm_bindgen(getter)]
        pub fn column(&self) -> usize {
            self.inner.column()
        }
        #[wasm_bindgen(getter)]
        pub fn category(&self) -> SerdeJsonCategory {
            match self.inner.classify() {
                serde_json::error::Category::Io => SerdeJsonCategory::Io,
                serde_json::error::Category::Syntax => SerdeJsonCategory::Syntax,
                serde_json::error::Category::Data => SerdeJsonCategory::Data,
                serde_json::error::Category::Eof => SerdeJsonCategory::Eof,
            }
        }
        #[wasm_bindgen(getter,js_name=ioErrorKind)]
        pub fn io_error_kind(&self) -> Option<crate::external_errors::IoErrorKind> {
            self.inner
                .io_error_kind()
                .map(crate::external_errors::IoErrorKind)
        }
    },
    source = crate::external_errors::StandardCause,
    cause = json_cause
);
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "SerdeJsonError")]
    #[derive(Clone)]
    pub type JsonCause;
    #[wasm_bindgen(
        typescript_type = "FingerprintError | FingerprintJsonError | SmartsWriteError | AtomPairError | FingerprintWorkerError | MatrixError | ValenceError | PeriodicTableError | LegacyStereoError | MoleculePropertyError | PropertyStringError | SmartsParseError | QueryMatchContextError | SubstructMatchError"
    )]
    #[derive(Clone)]
    pub type MorganCause;
    #[wasm_bindgen(
        typescript_type = "FingerprintError | AtomCodeError | MatrixError | TopologyValidationError | MorganError"
    )]
    #[derive(Clone)]
    pub type AtomPairCause;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    #[derive(Clone)]
    pub type PreparationCause;
    #[wasm_bindgen(typescript_type = "FingerprintPreparationError | MorganError")]
    #[derive(Clone)]
    pub type MorganReadCause;
}
macro_rules! child {
    ($owner:expr,$variant:path,$binding:path) => {
        <$binding>::from_owner($owner.child(|e| match e {
            $variant(v) => v,
            _ => unreachable!("immutable native error variant"),
        }))
        .unchecked_into()
    };
}
pub(crate) use child;
fn fingerprint_json_source(owner: &ReadonlyValue<ck::FingerprintJsonError>) -> Option<JsonCause> {
    match owner.as_ref() {
        ck::FingerprintJsonError::Parse(_) => Some(child!(
            owner,
            ck::FingerprintJsonError::Parse,
            SerdeJsonError
        )),
        _ => None,
    }
}
enum_error!(FingerprintJsonError,FingerprintJsonErrorKind,ck::FingerprintJsonError,{
    Parse=>ck::FingerprintJsonError::Parse(_),Invalid=>ck::FingerprintJsonError::Invalid(_),UnsupportedComponent=>ck::FingerprintJsonError::UnsupportedComponent{..}
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<JsonCause>{fingerprint_json_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn invalid(&self)->Option<String>{match self.inner.as_ref(){ck::FingerprintJsonError::Invalid(v)=>Some(v.clone()),_=>None}}
    #[wasm_bindgen(getter)] pub fn component(&self)->Option<String>{match self.inner.as_ref(){ck::FingerprintJsonError::UnsupportedComponent{component,..}=>Some((*component).into()),_=>None}}
    #[wasm_bindgen(getter,js_name=sourceType)] pub fn source_type(&self)->Option<String>{match self.inner.as_ref(){ck::FingerprintJsonError::UnsupportedComponent{source_type,..}=>Some(source_type.clone()),_=>None}}
},source=JsonCause,cause=fingerprint_json_source);
fn morgan_source(owner: &ReadonlyValue<ck::MorganError>) -> Option<MorganCause> {
    use ck::MorganError as E;
    Some(match owner.as_ref() {
        E::Fingerprint(_) => child!(
            owner,
            E::Fingerprint,
            crate::fingerprint_errors::FingerprintError
        ),
        E::Json(_) => child!(owner, E::Json, FingerprintJsonError),
        E::SmartsWrite(_) => child!(
            owner,
            E::SmartsWrite,
            crate::smarts_errors::SmartsWriteError
        ),
        E::AtomPair(_) => AtomPairError::from_owner(owner.child(|e| match e {
            E::AtomPair(v) => v.as_ref(),
            _ => unreachable!("immutable variant"),
        }))
        .unchecked_into(),
        E::Worker(_) => child!(
            owner,
            E::Worker,
            crate::worker_errors::FingerprintWorkerError
        ),
        E::StatePoisoned => return None,
        E::Matrix(_) => child!(owner, E::Matrix, crate::matrix_errors::MatrixError),
        E::Valence(_) => child!(owner, E::Valence, crate::valence_errors::ValenceError),
        E::PeriodicTable(_) => child!(
            owner,
            E::PeriodicTable,
            crate::periodic_errors::PeriodicTableError
        ),
        E::LegacyStereo(_) => child!(
            owner,
            E::LegacyStereo,
            crate::legacy_stereo_errors::LegacyStereoError
        ),
        E::MoleculeProperty(_) => child!(
            owner,
            E::MoleculeProperty,
            crate::model_errors::MoleculePropertyError
        ),
        E::PropertyString(_) => child!(
            owner,
            E::PropertyString,
            crate::property_errors::PropertyStringError
        ),
        E::SmartsParse(_) => child!(
            owner,
            E::SmartsParse,
            crate::smarts_parse_errors::SmartsParseError
        ),
        E::QueryMatchContext(_) => child!(
            owner,
            E::QueryMatchContext,
            crate::match_errors::QueryMatchContextError
        ),
        E::SubstructMatch(_) => child!(
            owner,
            E::SubstructMatch,
            crate::match_errors::SubstructMatchError
        ),
    })
}
enum_error!(MorganError,MorganErrorKind,ck::MorganError,{
    Fingerprint=>ck::MorganError::Fingerprint(_),Json=>ck::MorganError::Json(_),SmartsWrite=>ck::MorganError::SmartsWrite(_),AtomPair=>ck::MorganError::AtomPair(_),Worker=>ck::MorganError::Worker(_),StatePoisoned=>ck::MorganError::StatePoisoned,Matrix=>ck::MorganError::Matrix(_),Valence=>ck::MorganError::Valence(_),PeriodicTable=>ck::MorganError::PeriodicTable(_),LegacyStereo=>ck::MorganError::LegacyStereo(_),MoleculeProperty=>ck::MorganError::MoleculeProperty(_),PropertyString=>ck::MorganError::PropertyString(_),SmartsParse=>ck::MorganError::SmartsParse(_),QueryMatchContext=>ck::MorganError::QueryMatchContext(_),SubstructMatch=>ck::MorganError::SubstructMatch(_)
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<MorganCause>{morgan_source(&self.inner)}
},source=MorganCause,cause=morgan_source);
fn atom_pair_source(owner: &ReadonlyValue<ck::AtomPairError>) -> Option<AtomPairCause> {
    use ck::AtomPairError as E;
    Some(match owner.as_ref() {
        E::Fingerprint(_) => child!(
            owner,
            E::Fingerprint,
            crate::fingerprint_errors::FingerprintError
        ),
        E::AtomCode(_) => child!(owner, E::AtomCode, crate::fingerprint_errors::AtomCodeError),
        E::Matrix(_) => child!(owner, E::Matrix, crate::matrix_errors::MatrixError),
        E::Topology(_) => child!(
            owner,
            E::Topology,
            crate::model_errors::TopologyValidationError
        ),
        E::Preparation(_) => child!(owner, E::Preparation, MorganError),
        E::AtomInvariantLength { .. } => return None,
    })
}
enum_error!(AtomPairError,AtomPairErrorKind,ck::AtomPairError,{
    Fingerprint=>ck::AtomPairError::Fingerprint(_),AtomCode=>ck::AtomPairError::AtomCode(_),Matrix=>ck::AtomPairError::Matrix(_),Topology=>ck::AtomPairError::Topology(_),Preparation=>ck::AtomPairError::Preparation(_),AtomInvariantLength=>ck::AtomPairError::AtomInvariantLength{..}
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<AtomPairCause>{atom_pair_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn length(&self)->Option<usize>{match self.inner.as_ref(){ck::AtomPairError::AtomInvariantLength{length,..}=>Some(*length),_=>None}}
    #[wasm_bindgen(getter,js_name=atomIndex)] pub fn atom_index(&self)->Option<usize>{match self.inner.as_ref(){ck::AtomPairError::AtomInvariantLength{atom_index,..}=>Some(*atom_index),_=>None}}
},source=AtomPairCause,cause=atom_pair_source);
fn preparation_source(
    owner: &ReadonlyValue<ck::FingerprintPreparationError>,
) -> Option<PreparationCause> {
    match owner.as_ref() {
        ck::FingerprintPreparationError::MissingPreparedValence => None,
        ck::FingerprintPreparationError::RingPreparation(_) => Some(child!(
            owner,
            ck::FingerprintPreparationError::RingPreparation,
            crate::ring_errors::RingFindingError
        )),
    }
}
enum_error!(FingerprintPreparationError,FingerprintPreparationErrorKind,ck::FingerprintPreparationError,{
    MissingPreparedValence=>ck::FingerprintPreparationError::MissingPreparedValence,RingPreparation=>ck::FingerprintPreparationError::RingPreparation(_)
},{#[wasm_bindgen(getter)] pub fn source(&self)->Option<PreparationCause>{preparation_source(&self.inner)}},source=PreparationCause,cause=preparation_source);
fn read_source(owner: &ReadonlyValue<ck::MorganReadError>) -> Option<MorganReadCause> {
    Some(match owner.as_ref() {
        ck::MorganReadError::Preparation(_) => child!(
            owner,
            ck::MorganReadError::Preparation,
            FingerprintPreparationError
        ),
        ck::MorganReadError::Generator(_) => {
            child!(owner, ck::MorganReadError::Generator, MorganError)
        }
    })
}
enum_error!(MorganReadError,MorganReadErrorKind,ck::MorganReadError,{
    Preparation=>ck::MorganReadError::Preparation(_),Generator=>ck::MorganReadError::Generator(_)
},{#[wasm_bindgen(getter)] pub fn source(&self)->Option<MorganReadCause>{read_source(&self.inner)}},source=MorganReadCause,cause=read_source);
