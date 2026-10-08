//! Typed errors from canonical fingerprint values, without changing categories.
use crate::alignment_values::set;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct FingerprintError {
    inner: ck::FingerprintError,
}
#[wasm_bindgen]
impl FingerprintError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "fingerprints".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "never | null")]
    pub fn cause(&self) -> JsValue {
        JsValue::NULL
    }
}
fn kind(source: &ck::FingerprintError) -> &'static str {
    use ck::FingerprintError as E;
    match source {
        E::EmptyFingerprint => "EmptyFingerprint",
        E::Unsupported => "Unsupported",
        E::SparseIndexOutOfRange { .. } => "SparseIndexOutOfRange",
        E::BitLengthMismatch { .. } => "BitLengthMismatch",
        E::InvalidFoldFactor { .. } => "InvalidFoldFactor",
        E::RangeError { .. } => "RangeError",
        E::UndefinedArithmetic { .. } => "UndefinedArithmetic",
        E::PreconditionViolation { .. } => "PreconditionViolation",
        E::InvalidArguments { .. } => "InvalidArguments",
    }
}
pub(crate) fn fingerprint_error(source: &ck::FingerprintError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: canonical_values.rs::fingerprint_pyerr:
    //     let error = annotate(py, FingerprintError::new_err(source.to_string()), "fingerprints", kind, &source);
    // Boundary copies the small public error; no chemistry or arithmetic is reimplemented.
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("FingerprintError");
    let error: JsValue = error.into();
    set(&error, "domain", "fingerprints".into())?;
    set(&error, "kind", kind(source).into())?;
    set(&error, "detail", FingerprintError { inner: *source }.into())?;
    use ck::FingerprintError as E;
    match *source {
        E::SparseIndexOutOfRange { index, size } => {
            set(&error, "index", index.into())?;
            set(&error, "size", size.into())?;
        }
        E::BitLengthMismatch { left, right } => {
            set(&error, "left", left.into())?;
            set(&error, "right", right.into())?;
        }
        E::InvalidFoldFactor { factor, n_bits } => {
            set(&error, "factor", factor.into())?;
            set(&error, "nBits", n_bits.into())?;
        }
        E::RangeError { value } => set(&error, "value", value.into())?,
        E::UndefinedArithmetic { site } => set(&error, "site", site.into())?,
        E::PreconditionViolation { what } => set(&error, "what", what.into())?,
        E::InvalidArguments { reason } => set(&error, "reason", reason.into())?,
        E::EmptyFingerprint | E::Unsupported => (),
    }
    Ok(error)
}
pub(crate) fn throw_fingerprint_error(source: ck::FingerprintError) -> JsValue {
    fingerprint_error(&source).unwrap_or_else(|e| e)
}
