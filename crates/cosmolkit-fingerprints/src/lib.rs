//! Detached molecular fingerprint boundaries.
//!
//! Value layer: source-backed ports of RDKit `Code/DataStructs`
//! `ExplicitBitVect`, `SparseBitVect`, `SparseIntVect`, `BitOps` similarity
//! and folding, and the `RDGeneral/hash` integral closure used by
//! fingerprints. Generator entrypoints remain explicit unsupported
//! capability boundaries.

pub mod folding;
pub mod hash;
pub mod similarity;
mod sparse_bits;
mod sparse_counts;
mod values;

use std::fmt;

use cosmolkit_model::TopologyBlock;

pub use sparse_bits::SparseBitFingerprint;
pub use sparse_counts::{SparseCountFingerprint, SparseCountFingerprint32};
pub use values::Fingerprint;

/// Structured fingerprint value/generator errors.
///
/// Error categories mirror the pinned RDKit exception categories
/// (`IndexErrorException`, `ValueErrorException`, `RangeErrorException`)
/// plus the existing fail-closed `Unsupported` generator boundary.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FingerprintError {
    Unsupported,
    /// Source `IndexErrorException`: an index is outside the vector length.
    SparseIndexOutOfRange {
        index: u64,
        size: u64,
    },
    /// Source `ValueErrorException("BitVects must be same length")` /
    /// `("SparseIntVect size mismatch")`. Fields are u64 so count-vector
    /// lengths (SparseIntVect<IndexType> with 64-bit IndexType) are
    /// reported without lossy narrowing.
    BitLengthMismatch {
        left: u64,
        right: u64,
    },
    /// Source `ValueErrorException("invalid fold factor")`.
    InvalidFoldFactor {
        factor: u32,
        n_bits: u32,
    },
    /// Source `RangeErrorException` from `RANGE_CHECK(0, a, 1)` /
    /// `RANGE_CHECK(0, b, 1)` in `TverskySimilarity`.
    RangeError {
        value: f64,
    },
    /// Structured rejection at an actually executed signed-i32 operation
    /// whose pinned C++ source is undefined (overflow, `abs(MIN)`,
    /// executed `/0`, `MIN / -1`, right-only negation of MIN). Per the
    /// Supervisor resolution of 2026-09-24 this is a COSMolKit safety
    /// boundary — not an RDKit exception and not a claim of all-input
    /// parity. Defined source arithmetic remains unchanged.
    UndefinedArithmetic {
        site: &'static str,
    },
    /// Source `Invar::Invariant` precondition violation (for example the
    /// bitmap fast path of `NumOnBitsInCommon<ExplicitBitVect>` on a
    /// zero-length vector, `PRECONDITION(afp, "no afp")`,
    /// BitOps.cpp:954). Reproduced boundary, oracle-verified.
    PreconditionViolation {
        what: &'static str,
    },
    /// Generic invalid-arguments surface retaining source-defined rejection
    /// text and explicitly documented safe rejections of undefined operations.
    InvalidArguments {
        reason: &'static str,
    },
}

// No `impl Eq`: `RangeError` carries an `f64`, whose derived `PartialEq`
// does not admit a total-equivalence promise (F54 correction of an
// earlier unsound manual `impl Eq`). Consumer audit: no current consumer
// inside this crate (or in the workspace, per inventory.md — zero Rust
// consumers) requires `FingerprintError: Eq`; any future protected
// consumer needing Eq is a supervisor decision, not a local change.

impl fmt::Display for FingerprintError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unsupported => {
                write!(f, "unsupported fingerprint operation")
            }
            Self::SparseIndexOutOfRange { index, size } => {
                write!(
                    f,
                    "fingerprint index {index} is outside vector length {size}"
                )
            }
            Self::BitLengthMismatch { left, right } => {
                write!(f, "fingerprint bit length mismatch: {left} != {right}")
            }
            Self::InvalidFoldFactor { factor, n_bits } => {
                write!(
                    f,
                    "invalid fold factor {factor} for fingerprint of {n_bits} bits"
                )
            }
            Self::RangeError { value } => {
                write!(f, "similarity parameter outside [0,1]: {value}")
            }
            Self::InvalidArguments { reason } => {
                write!(f, "invalid fingerprint arguments: {reason}")
            }
            Self::UndefinedArithmetic { site } => {
                write!(
                    f,
                    "undefined signed arithmetic at {site} (COSMolKit safety \
                     boundary on a C++-undefined operation; not an RDKit \
                     exception)"
                )
            }
            Self::PreconditionViolation { what } => {
                write!(f, "pre-condition violation: {what}")
            }
        }
    }
}

impl std::error::Error for FingerprintError {}

pub fn morgan(topology: &TopologyBlock) -> Result<Fingerprint, FingerprintError> {
    let _ = topology;
    Err(FingerprintError::Unsupported)
}

pub fn pattern(topology: &TopologyBlock) -> Result<Fingerprint, FingerprintError> {
    morgan(topology)
}
