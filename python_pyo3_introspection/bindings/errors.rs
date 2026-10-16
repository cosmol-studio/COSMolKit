use crate::error_framework::{error_projection, variant_pattern};
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;
error_projection!(SmilesError, SmilesErrorKind, ck::SmilesError, "smiles", {
    Parse(_), Hydrogen(_), Sanitize(_), Stereo(_), Construction(_),
});

error_projection!(MaccsFingerprintError, MaccsFingerprintErrorKind, ck::MaccsFingerprintError, "fingerprint", {
    UnsupportedOption { option: &'static str, reason: &'static str },
    MissingPattern { bit: usize },
    Topology(_), Rings(_), Valence(_), Paths(_), Smarts(_), QueryCompile(_),
    QueryContext(_), Match(_), Value(_),
});

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
