//! Worker failures remain native diagnostics, not claims of browser threading.
use crate::errors::{enum_error, exception};
use crate::storage::ReadonlyValue;
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "IoError")]
    #[derive(Clone)]
    pub type ThreadCountCause;
    #[wasm_bindgen(typescript_type = "ThreadCountError | IoError")]
    #[derive(Clone)]
    pub type FingerprintWorkerCause;
}
fn count_source(owner: &ReadonlyValue<ck::ThreadCountError>) -> Option<ThreadCountCause> {
    match owner.as_ref() {
        ck::ThreadCountError::ObservationIo(_) => Some(
            crate::external_errors::IoError::from_owner(owner.child(|e| match e {
                ck::ThreadCountError::ObservationIo(v) => v,
                _ => unreachable!("immutable selected IO variant"),
            }))
            .unchecked_into(),
        ),
        _ => None,
    }
}
enum_error!(ThreadCountError,ThreadCountErrorKind,ck::ThreadCountError,{
    UndefinedSignedNegation=>ck::ThreadCountError::UndefinedSignedNegation,ObservationIo=>ck::ThreadCountError::ObservationIo(_),InvalidCpuList=>ck::ThreadCountError::InvalidCpuList{..},AutomaticObservationUnavailable=>ck::ThreadCountError::AutomaticObservationUnavailable{..},
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<ThreadCountCause>{count_source(&self.inner)}
    #[wasm_bindgen(getter)] pub fn segment(&self)->Option<usize>{match self.inner.as_ref(){ck::ThreadCountError::InvalidCpuList{segment,..}=>Some(*segment),_=>None}}
    #[wasm_bindgen(getter)] pub fn reason(&self)->Option<String>{match self.inner.as_ref(){ck::ThreadCountError::InvalidCpuList{reason,..}=>Some((*reason).into()),_=>None}}
    #[wasm_bindgen(getter)] pub fn platform(&self)->Option<String>{match self.inner.as_ref(){ck::ThreadCountError::AutomaticObservationUnavailable{platform}=>Some((*platform).into()),_=>None}}
},source=ThreadCountCause,cause=count_source);
#[cfg(feature = "cap-fingerprints")]
fn worker_source(
    owner: &ReadonlyValue<ck::FingerprintWorkerError>,
) -> Option<FingerprintWorkerCause> {
    use ck::FingerprintWorkerError as E;
    Some(
        match owner.as_ref() {
            E::ThreadCount(_) => ThreadCountError::from_owner(owner.child(|e| match e {
                E::ThreadCount(v) => v,
                _ => unreachable!("immutable selected count variant"),
            })),
            E::ThreadSpawn(_) => {
                crate::external_errors::IoError::from_owner(owner.child(|e| match e {
                    E::ThreadSpawn(v) => v,
                    _ => unreachable!("immutable selected spawn variant"),
                }))
            }
            _ => return None,
        }
        .unchecked_into(),
    )
}
#[cfg(feature = "cap-fingerprints")]
enum_error!(FingerprintWorkerError,FingerprintWorkerErrorKind,ck::FingerprintWorkerError,{
    ThreadCount=>ck::FingerprintWorkerError::ThreadCount(_),ThreadSpawn=>ck::FingerprintWorkerError::ThreadSpawn(_),Panic=>ck::FingerprintWorkerError::Panic,Protocol=>ck::FingerprintWorkerError::Protocol,
},{
    #[wasm_bindgen(getter)] pub fn source(&self)->Option<FingerprintWorkerCause>{worker_source(&self.inner)}
},source=FingerprintWorkerCause,cause=worker_source);
