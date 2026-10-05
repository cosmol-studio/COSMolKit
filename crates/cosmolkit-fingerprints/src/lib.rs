//! Detached molecular fingerprint boundaries.
//!
//! Value layer: source-backed ports of RDKit `Code/DataStructs`
//! `ExplicitBitVect`, `SparseBitVect`, `SparseIntVect`, `BitOps` similarity
//! and folding, and the `RDGeneral/hash` integral closure used by
//! fingerprints. Narrow Morgan generator calls operate on explicitly
//! prepared detached molecule state.
//!
//! The detached Morgan surface exposes canonical call/configuration values
//! and the four completed result functions while keeping implementation
//! modules private. These function-pointer assignments compile against their
//! exact exported signatures:
//!
//! ```rust
//! use cosmolkit_fingerprints::{
//!     FingerprintAdditionalOutput, Fingerprint, MorganAtomInvariants, MorganCall, MorganError,
//!     MorganParams, MorganPreparedInput, SparseBitFingerprint, SparseCountFingerprint,
//!     SparseCountFingerprint32, morgan_bits, morgan_count, morgan_sparse_bits,
//!     morgan_sparse_count,
//! };
//!
//! fn main() {
//!     let _: fn(
//!         &MorganPreparedInput<'_>,
//!         &MorganParams,
//!         &MorganCall<'_>,
//!         MorganAtomInvariants<'_>,
//!         Option<&mut FingerprintAdditionalOutput>,
//!     ) -> Result<SparseCountFingerprint, MorganError> = morgan_sparse_count;
//!     let _: fn(
//!         &MorganPreparedInput<'_>,
//!         &MorganParams,
//!         &MorganCall<'_>,
//!         MorganAtomInvariants<'_>,
//!         Option<&mut FingerprintAdditionalOutput>,
//!     ) -> Result<SparseBitFingerprint, MorganError> = morgan_sparse_bits;
//!     let _: fn(
//!         &MorganPreparedInput<'_>,
//!         &MorganParams,
//!         &MorganCall<'_>,
//!         MorganAtomInvariants<'_>,
//!         Option<&mut FingerprintAdditionalOutput>,
//!     ) -> Result<SparseCountFingerprint32, MorganError> = morgan_count;
//!     let _: fn(
//!         &MorganPreparedInput<'_>,
//!         &MorganParams,
//!         &MorganCall<'_>,
//!         MorganAtomInvariants<'_>,
//!         Option<&mut FingerprintAdditionalOutput>,
//!     ) -> Result<Fingerprint, MorganError> = morgan_bits;
//! }
//! ```
//!
//! Implementation modules are not part of that boundary:
//!
//! ```compile_fail
//! use cosmolkit_fingerprints::morgan::MorganGenerator;
//! ```
//!
//! ```compile_fail
//! use cosmolkit_fingerprints::prepared;
//! ```
//!
//! ```compile_fail
//! use cosmolkit_fingerprints::generator;
//! ```
//!
//! ```compile_fail
//! use cosmolkit_fingerprints::rng;
//! ```

mod additional_output;
mod atom_code;
mod atom_pair;
pub mod folding;
mod generator;
pub mod hash;
mod invariants;
mod layered;
mod molecule_hash;
mod morgan;
mod packed_codes;
mod prepared;
mod rng;
pub mod similarity;
mod sparse_bits;
mod sparse_counts;
mod topological_torsion;
mod values;

use std::fmt;

use cosmolkit_core::{
    LegacyStereoError, MatrixError, PeriodicTableError, PropertyStringError, ValenceError,
};
use cosmolkit_model::{MoleculePropertyError, TopologyBlock};

pub use additional_output::FingerprintAdditionalOutput;
pub use atom_code::{AtomCodeAssignment, AtomCodeError, AtomCodeInput, AtomCodeOptions, atom_code};
pub use atom_pair::{
    AtomPairAtomInvariantsGenerator, AtomPairCall, AtomPairError, AtomPairParams,
    AtomPairPreparedInput, atom_pair_bits, atom_pair_count, atom_pair_sparse_bits,
    atom_pair_sparse_count,
};
pub use cosmolkit_search::QueryGraph;
pub use layered::{
    LAYERED_FINGERPRINT_MAX_LAYERS, LAYERED_FINGERPRINT_SUBSTRUCTURE_LAYERS,
    LAYERED_FINGERPRINT_VERSION, LayeredFingerprintError, LayeredFingerprintLayers,
    LayeredFingerprintParams, LayeredFingerprintResult, layered_fingerprint,
    layered_fingerprint_with_output, layered_query_fingerprint,
    layered_query_fingerprint_with_output,
};
pub use morgan::{
    MorganAtomInvariants, MorganCall, MorganParams, morgan_bits, morgan_count, morgan_sparse_bits,
    morgan_sparse_count,
};
pub use packed_codes::{atom_pair_code, topological_torsion_code, topological_torsion_hash};
pub use prepared::MorganPreparedInput;
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

/// Error causes shared by the detached Morgan generator implementation.
///
/// `FingerprintError` remains the small Copy value error used by the existing
/// fingerprint-value API; chemistry and table failures retain their concrete
/// source types here instead of being stringified.
#[derive(Debug)]
pub enum MorganError {
    Fingerprint(FingerprintError),
    Json(FingerprintJsonError),
    SmartsWrite(cosmolkit_search::SmartsWriteError),
    AtomPair(Box<AtomPairError>),
    Worker(FingerprintWorkerError),
    StatePoisoned,
    Matrix(MatrixError),
    Valence(ValenceError),
    PeriodicTable(PeriodicTableError),
    LegacyStereo(LegacyStereoError),
    MoleculeProperty(MoleculePropertyError),
    PropertyString(PropertyStringError),
    SmartsParse(cosmolkit_search::SmartsParseError),
    QueryMatchContext(cosmolkit_search::QueryMatchContextError),
    SubstructMatch(cosmolkit_search::SubstructMatchError),
}

impl fmt::Display for MorganError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Fingerprint(source) => write!(f, "Morgan fingerprint error: {source}"),
            Self::Json(source) => source.fmt(f),
            Self::SmartsWrite(source) => source.fmt(f),
            Self::AtomPair(source) => source.fmt(f),
            Self::Worker(source) => source.fmt(f),
            Self::StatePoisoned => f.write_str("Morgan generator state lock was poisoned"),
            Self::Matrix(source) => write!(f, "Morgan matrix error: {source}"),
            Self::Valence(source) => write!(f, "Morgan valence preparation error: {source}"),
            Self::PeriodicTable(source) => {
                write!(f, "Morgan periodic-table lookup error: {source}")
            }
            Self::LegacyStereo(source) => {
                write!(f, "Morgan legacy stereo preparation error: {source}")
            }
            Self::MoleculeProperty(source) => {
                write!(f, "Morgan molecule-property error: {source}")
            }
            Self::PropertyString(source) => {
                write!(f, "Morgan property string conversion error: {source}")
            }
            Self::SmartsParse(source) => write!(f, "Morgan feature SMARTS parse error: {source}"),
            Self::QueryMatchContext(source) => {
                write!(f, "Morgan query-context validation error: {source}")
            }
            Self::SubstructMatch(source) => write!(f, "Morgan feature matching error: {source}"),
        }
    }
}

impl std::error::Error for MorganError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Json(source) => Some(source),
            Self::SmartsWrite(source) => Some(source),
            Self::AtomPair(source) => Some(source.as_ref()),
            Self::Worker(source) => Some(source),
            Self::StatePoisoned => None,
            Self::Fingerprint(source) => Some(source),
            Self::Matrix(source) => Some(source),
            Self::Valence(source) => Some(source),
            Self::PeriodicTable(source) => Some(source),
            Self::LegacyStereo(source) => Some(source),
            Self::MoleculeProperty(source) => Some(source),
            Self::PropertyString(source) => Some(source),
            Self::SmartsParse(source) => Some(source),
            Self::QueryMatchContext(source) => Some(source),
            Self::SubstructMatch(source) => Some(source),
        }
    }
}

impl From<FingerprintError> for MorganError {
    fn from(source: FingerprintError) -> Self {
        Self::Fingerprint(source)
    }
}

impl From<MatrixError> for MorganError {
    fn from(source: MatrixError) -> Self {
        Self::Matrix(source)
    }
}

impl From<ValenceError> for MorganError {
    fn from(source: ValenceError) -> Self {
        Self::Valence(source)
    }
}

impl From<PeriodicTableError> for MorganError {
    fn from(source: PeriodicTableError) -> Self {
        Self::PeriodicTable(source)
    }
}

impl From<LegacyStereoError> for MorganError {
    fn from(source: LegacyStereoError) -> Self {
        Self::LegacyStereo(source)
    }
}

impl From<MoleculePropertyError> for MorganError {
    fn from(source: MoleculePropertyError) -> Self {
        Self::MoleculeProperty(source)
    }
}

impl From<PropertyStringError> for MorganError {
    fn from(source: PropertyStringError) -> Self {
        Self::PropertyString(source)
    }
}

impl From<cosmolkit_search::SmartsParseError> for MorganError {
    fn from(source: cosmolkit_search::SmartsParseError) -> Self {
        Self::SmartsParse(source)
    }
}

impl From<cosmolkit_search::QueryMatchContextError> for MorganError {
    fn from(source: cosmolkit_search::QueryMatchContextError) -> Self {
        Self::QueryMatchContext(source)
    }
}

impl From<cosmolkit_search::SubstructMatchError> for MorganError {
    fn from(source: cosmolkit_search::SubstructMatchError) -> Self {
        Self::SubstructMatch(source)
    }
}

pub fn morgan(topology: &TopologyBlock) -> Result<Fingerprint, FingerprintError> {
    let _ = topology;
    Err(FingerprintError::Unsupported)
}

pub fn pattern(topology: &TopologyBlock) -> Result<Fingerprint, FingerprintError> {
    morgan(topology)
}

pub use topological_torsion::{
    TopologicalTorsionCall, TopologicalTorsionError, TopologicalTorsionParams,
    topological_torsion_bits, topological_torsion_count, topological_torsion_sparse_bits,
    topological_torsion_sparse_count,
};

#[cfg(test)]
mod test_support;

mod metadata;
pub use metadata::FingerprintJsonError;

pub use topological_torsion::{
    LegacyTopologicalTorsionParams, legacy_topological_torsion_bits,
    legacy_topological_torsion_count, legacy_topological_torsion_sparse_count,
};

pub use topological_torsion::{TopologicalTorsionGenerator, TopologicalTorsionSettings};

pub use topological_torsion::topological_torsion_ids;

mod atom_code_explanation;
pub use atom_code_explanation::{AtomCodeExplanation, AtomCodeExplanationError};
pub use molecule_hash::{MoleculeHashError, molecule_hash, molecule_hash_with_ranks};

mod argument_metadata;

mod fingerprint_bulk;
pub use fingerprint_bulk::FingerprintWorkerError;
pub use morgan::{MorganAtomProvider, MorganBondProvider, MorganOperator, MorganSettings};
impl From<FingerprintJsonError> for MorganError {
    fn from(e: FingerprintJsonError) -> Self {
        Self::Json(e)
    }
}
impl From<cosmolkit_search::SmartsWriteError> for MorganError {
    fn from(e: cosmolkit_search::SmartsWriteError) -> Self {
        Self::SmartsWrite(e)
    }
}
impl From<AtomPairError> for MorganError {
    fn from(e: AtomPairError) -> Self {
        Self::AtomPair(Box::new(e))
    }
}
impl From<FingerprintWorkerError> for MorganError {
    fn from(e: FingerprintWorkerError) -> Self {
        Self::Worker(e)
    }
}

mod atom_pairs_parameters;
pub use atom_pairs_parameters::AtomPairsParameters;

mod topological_torsion_path;
pub use topological_torsion_path::{
    TopologicalTorsionPathScoreError, explain_path_score, topological_torsion_path_score,
};

mod maccs;
pub use maccs::{
    MaccsFingerprintError, MaccsFingerprintParams, maccs_fingerprint, maccs_fingerprint_raw,
};

mod pattern;
pub use pattern::{
    PATTERN_FINGERPRINT_VERSION, PatternFingerprintError, PatternFingerprintParams,
    pattern_fingerprint, pattern_query_fingerprint,
};
