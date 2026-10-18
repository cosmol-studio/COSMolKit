use crate::error_framework::error_projection;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;
error_projection!(AtomCodeExplanationError, AtomCodeExplanationErrorKind, ck::AtomCodeExplanationError, "fingerprint", {
    UnknownChirality { code: u8 },
});
error_projection!(Coordinate3DReadError, Coordinate3DReadErrorKind, ck::Coordinate3DReadError, "coordinates", {
    ConformerNotFound { conformer_id: usize, count: usize },
});
pub use crate::operation_errors::{SmilesError, SmilesErrorKind};

pub use crate::maccs_errors::{MaccsFingerprintError, MaccsFingerprintErrorKind};

error_projection!(FingerprintError, FingerprintErrorKind, ck::FingerprintError, "fingerprint", {
    EmptyFingerprint,
    Unsupported,
    SparseIndexOutOfRange { index: u64, size: u64 },
    BitLengthMismatch { left: u64, right: u64 },
    InvalidFoldFactor { factor: u32, n_bits: u32 },
    RangeError { value: f64 },
    UndefinedArithmetic { site: &'static str },
    PreconditionViolation { what: &'static str },
    InvalidArguments { reason: &'static str },
});
