//! Standard errors retain their actual native owner and source chain.
use crate::errors::{error_projection, exception};
use crate::storage::ReadonlyValue;
use std::error::Error;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "IoError | Utf8Error | FromUtf8Error | ParseFloatError | ExternalError"
    )]
    #[derive(Clone)]
    pub type StandardCause;
    #[wasm_bindgen(typescript_type = "Utf8Error")]
    #[derive(Clone)]
    pub type Utf8Cause;
}

error_projection!(Utf8Error, std::str::Utf8Error, {
    #[wasm_bindgen(getter,js_name=validUpTo)]
    pub fn valid_up_to(&self) -> usize {
        self.inner.valid_up_to()
    }
    #[wasm_bindgen(getter,js_name=errorLen)]
    pub fn error_len(&self) -> Option<usize> {
        self.inner.error_len()
    }
});
error_projection!(FromUtf8Error, std::string::FromUtf8Error, {
    #[wasm_bindgen(getter)]
    pub fn bytes(&self) -> Vec<u8> {
        self.inner.as_bytes().to_vec()
    }
    #[wasm_bindgen(getter,js_name=utf8Error)]
    pub fn utf8_error(&self) -> Utf8Cause {
        Utf8Error::from_rust(self.inner.utf8_error()).unchecked_into()
    }
});

#[wasm_bindgen]
pub struct IoErrorKind(pub(crate) std::io::ErrorKind);
#[wasm_bindgen]
impl IoErrorKind {
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        format!("{:?}", self.0)
    }
    #[wasm_bindgen(js_name=toString)]
    pub fn display(&self) -> String {
        self.0.to_string()
    }
}

fn io_cause(owner: &ReadonlyValue<std::io::Error>) -> Option<StandardCause> {
    owner.source().map(|_| {
        standard_source(owner.child(|value| value.source().expect("existing immutable IO source")))
            .unchecked_into()
    })
}
error_projection!(
    IoError,
    std::io::Error,
    {
        #[wasm_bindgen(getter)]
        pub fn kind(&self) -> IoErrorKind {
            IoErrorKind(self.inner.kind())
        }
        #[wasm_bindgen(getter,js_name=rawOsError)]
        pub fn raw_os_error(&self) -> Option<i32> {
            self.inner.raw_os_error()
        }
    },
    source = StandardCause,
    cause = io_cause
);

// Error::source is open-ended for external libraries. Known standard values are
// projected exactly; an unknown external type keeps its original Display and
// source chain. This is never used in place of a known CK error projection.
pub(crate) fn standard_source(owner: ReadonlyValue<dyn Error>) -> JsValue {
    if owner.is::<std::io::Error>() {
        return IoError::from_owner(owner.child(|e| e.downcast_ref().expect("checked IO type")));
    }
    if owner.is::<std::str::Utf8Error>() {
        return Utf8Error::from_owner(
            owner.child(|e| e.downcast_ref().expect("checked UTF-8 type")),
        );
    }
    if owner.is::<std::string::FromUtf8Error>() {
        return FromUtf8Error::from_owner(
            owner.child(|e| e.downcast_ref().expect("checked owned UTF-8 type")),
        );
    }
    if owner.is::<std::num::ParseFloatError>() {
        return crate::xyz_errors::ParseFloatError::from_owner(
            owner.child(|e| e.downcast_ref().expect("checked float-parser type")),
        );
    }
    ExternalError::from_owner(owner)
}

#[wasm_bindgen]
pub struct ExternalError {
    inner: ReadonlyValue<dyn Error>,
    cause: Option<StandardCause>,
}
impl ExternalError {
    fn from_owner(inner: ReadonlyValue<dyn Error>) -> JsValue {
        let cause = inner.source().map(|_| {
            standard_source(
                inner.child(|e| e.source().expect("existing immutable external source")),
            )
            .unchecked_into()
        });
        exception(Self { inner, cause }.into())
    }
}
#[wasm_bindgen]
impl ExternalError {
    #[wasm_bindgen(getter)]
    pub fn name(&self) -> String {
        "ExternalError".into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter)]
    pub fn cause(&self) -> Option<StandardCause> {
        self.cause.clone()
    }
    #[wasm_bindgen(js_name=toString)]
    pub fn display(&self) -> String {
        self.inner.to_string()
    }
}
