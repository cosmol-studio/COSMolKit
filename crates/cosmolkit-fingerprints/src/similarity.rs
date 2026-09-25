//! Scalar bit-vector similarity and counting operations: source-backed
//! port of the pinned RDKit `Code/DataStructs/BitOps.{h,cpp}` metric
//! family over the dense and sparse bit vectors.

use crate::{Fingerprint, FingerprintError, SparseBitFingerprint};

fn dense_same_length(a: &Fingerprint, b: &Fingerprint) -> Result<(), FingerprintError> {
    if a.n_bits() != b.n_bits() {
        return Err(FingerprintError::BitLengthMismatch {
            left: u64::from(a.n_bits()),
            right: u64::from(b.n_bits()),
        });
    }
    Ok(())
}

fn sparse_same_length(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<(), FingerprintError> {
    if a.n_bits() != b.n_bits() {
        return Err(FingerprintError::BitLengthMismatch {
            left: u64::from(a.n_bits()),
            right: u64::from(b.n_bits()),
        });
    }
    Ok(())
}

/// `NumOnBitsInCommon<ExplicitBitVect>` (BitOps.cpp:259-269).
///
/// The dense specialization computes `((bv1.dp_bits) & (bv2.dp_bits))
/// .count()`; the underlying boost operator raises on size mismatch,
/// mapped to the structured `BitLengthMismatch` error.
pub fn num_on_bits_in_common(a: &Fingerprint, b: &Fingerprint) -> Result<u32, FingerprintError> {
    // RDKit✔️✔️: int NumOnBitsInCommon(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2) {
    // RDKit✔️✔️:   // Don't try this at home, we (hope we) know what we're doing
    // RDKit✔️✔️:   const unsigned char *afp, *bfp;
    // RDKit✔️✔️:   unsigned int nBytes;
    // RDKit✔️✔️:   if (EBVToBitmap(bv1, afp, nBytes) && EBVToBitmap(bv2, bfp, nBytes)) {
    // RDKit✔️✔️:     unsigned int result = CalcBitmapNumBitsInCommon(afp, bfp, nBytes);
    // RDKit✔️✔️:     return (int)result;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<int>(((*bv1.dp_bits) & (*bv2.dp_bits)).count());
    // RDKit✔️✔️: }
    //
    // The pointer-cast bitmap fast path is an implementation detail of the
    // same AND-popcount result — except for zero-length vectors: the
    // empty `m_bits` storage yields a null bitmap pointer and
    // `CalcBitmapNumBitsInCommon`'s `PRECONDITION(afp, "no afp")` raises
    // (oracle-verified: every dense on-count metric errors on a
    // zero-length pair except Tanimoto, whose `total == 0` guard fires
    // first). Reproduced as a structured precondition error; non-empty
    // inputs take the counted-AND path with identical output.
    if a.n_bits() == 0 {
        return Err(FingerprintError::PreconditionViolation {
            what: "no afp (CalcBitmapNumBitsInCommon on zero-length vector)",
        });
    }
    Ok(a.and(b)?.num_on_bits())
}

/// `NumOnBitsInCommon<SparseBitVect>` — the generic template (BitOps.cpp
/// :222-225) delegates to `OnBitsInCommon`, whose length precondition
/// raises `ValueErrorException("BitVects must be same length")` before
/// any counting.
pub fn sparse_num_on_bits_in_common(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<u32, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: int NumOnBitsInCommon(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   return static_cast<int>(OnBitsInCommon(bv1, bv2).size());
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    Ok(a.bit_set().intersection(b.bit_set()).count() as u32)
}

/// `OnBitsInCommon` (BitOps.cpp:561-569): indices of bits on in both
/// vectors, in ascending order. Length precondition:
/// `ValueErrorException("BitVects must be same length")`.
pub fn on_bits_in_common(a: &Fingerprint, b: &Fingerprint) -> Result<Vec<u32>, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: IntVect OnBitsInCommon(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   IntVect res;
    // RDKit✔️✔️:   (bv1 & bv2).getOnBits(res);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    Ok(a.and(b)?.on_bits())
}

/// `OnBitsInCommon` over the sparse bit vectors; the returned indices are
/// the raw signed carrier keys in signed ascending order (F02), exactly
/// like the source's `getOnBits` on the intersected set.
pub fn sparse_on_bits_in_common(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<Vec<i32>, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: IntVect OnBitsInCommon(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   IntVect res;
    // RDKit✔️✔️:   (bv1 & bv2).getOnBits(res);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    Ok(a.intersection(b).on_bits())
}

/// `OffBitsInCommon` (BitOps.cpp:585-593): indices of bits off in both
/// vectors — `(~(bv1 | bv2)).getOnBits()` — with the length
/// precondition of the generic template.
pub fn off_bits_in_common(a: &Fingerprint, b: &Fingerprint) -> Result<Vec<u32>, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: IntVect OffBitsInCommon(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   IntVect res;
    // RDKit✔️✔️:   (~(bv1 | bv2)).getOnBits(res);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    Ok(a.or(b)?.not().on_bits())
}

/// `OffBitsInCommon` over the sparse bit vectors (raw signed keys, signed
/// order). The complement loop only visits indices below the length, so
/// retained out-of-length members of the union do not appear.
pub fn sparse_off_bits_in_common(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<Vec<i32>, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: IntVect OffBitsInCommon(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   IntVect res;
    // RDKit✔️✔️:   (~(bv1 | bv2)).getOnBits(res);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    Ok(a.union(b).negate().on_bits())
}

/// `NumBitsInCommon` (BitOps.cpp:510-522): number of bits (on and off)
/// that agree — `bv1.getNumBits() - (bv1 ^ bv2).getNumOnBits()`, with the
/// generic length precondition.
pub fn num_bits_in_common(a: &Fingerprint, b: &Fingerprint) -> Result<u32, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: int NumBitsInCommon(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return bv1.getNumBits() - (bv1 ^ bv2).getNumOnBits();
    // RDKit✔️✔️: }
    //
    // Dense specialization (BitOps.cpp:519-522):
    // RDKit✔️✔️: int NumBitsInCommon(const ExplicitBitVect &bv1, const ExplicitBitVect &bv2) {
    // RDKit✔️✔️:   return bv1.getNumBits() -
    // RDKit✔️✔️:          static_cast<int>(((*bv1.dp_bits) ^ (*bv2.dp_bits)).count());
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    Ok(a.n_bits() - a.xor(b)?.num_on_bits())
}

/// `NumBitsInCommon` over the sparse bit vectors.
pub fn sparse_num_bits_in_common(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<u32, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: int NumBitsInCommon(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return bv1.getNumBits() - (bv1 ^ bv2).getNumOnBits();
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    // The subtraction is C++ `unsigned int` arithmetic: with retained
    // out-of-length members the xor on-count can exceed the length, and
    // the defined source result wraps modulo 2^32 (reproduced, not an
    // error; pinned by regression).
    Ok(a.n_bits()
        .wrapping_sub(a.symmetric_difference(b).num_on_bits()))
}

/// `TanimotoSimilarity` (BitOps.cpp:286-297): with the source's exact
/// branch order — length precondition, `total == 0 -> 0.0` (both empty),
/// then `common / (total - common)`.
pub fn tanimoto(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double TanimotoSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int total = bv1.getNumOnBits() + bv2.getNumOnBits();
    // RDKit✔️✔️:   if (total == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int common = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   return (double)common / (double)(total - common);
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    let total = tanimoto_totals(a.num_on_bits(), b.num_on_bits());
    if total == 0 {
        return Ok(0.0);
    }
    let common = num_on_bits_in_common(a, b)?;
    Ok(f64::from(common) / f64::from(tanimoto_denominator(total, common)))
}

/// `TanimotoSimilarity` over the sparse bit vectors (same source body).
pub fn sparse_tanimoto(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double TanimotoSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int total = bv1.getNumOnBits() + bv2.getNumOnBits();
    // RDKit✔️✔️:   if (total == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int common = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   return (double)common / (double)(total - common);
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    let total = tanimoto_totals(a.num_on_bits(), b.num_on_bits());
    if total == 0 {
        return Ok(0.0);
    }
    let common = sparse_num_on_bits_in_common(a, b)?;
    Ok(f64::from(common) / f64::from(tanimoto_denominator(total, common)))
}

/// `TverskySimilarity` (BitOps.cpp:299-318) with the source's exact
/// branch order: `RANGE_CHECK(0,a,1)` then `RANGE_CHECK(0,b,1)` (raising
/// the source's Range Error invariant before the length precondition),
/// length check, `y == 0 || z == 0 -> 0.0`, then
/// `denom == 0.0 -> 1.0 else x / denom`.
pub fn tversky(
    a: &Fingerprint,
    b: &Fingerprint,
    alpha: f64,
    beta: f64,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double TverskySimilarity(const T1 &bv1, const T2 &bv2, double a, double b) {
    // RDKit✔️✔️:   RANGE_CHECK(0, a, 1);
    // RDKit✔️✔️:   RANGE_CHECK(0, b, 1);
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   auto y = bv1.getNumOnBits();
    // RDKit✔️✔️:   auto z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y == 0 || z == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double denom = a * y + b * z + (1 - a - b) * x;
    // RDKit✔️✔️:   if (denom == 0.0) {
    // RDKit✔️✔️:     return 1.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return x / denom;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    range_check(alpha)?;
    range_check(beta)?;
    dense_same_length(a, b)?;
    tversky_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
        alpha,
        beta,
    )
}

/// `TverskySimilarity` over the sparse bit vectors (same source body).
#[allow(clippy::similar_names)]
pub fn sparse_tversky(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
    alpha: f64,
    beta: f64,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double TverskySimilarity(const T1 &bv1, const T2 &bv2, double a, double b) {
    // RDKit✔️✔️:   RANGE_CHECK(0, a, 1);
    // RDKit✔️✔️:   RANGE_CHECK(0, b, 1);
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   auto y = bv1.getNumOnBits();
    // RDKit✔️✔️:   auto z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y == 0 || z == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double denom = a * y + b * z + (1 - a - b) * x;
    // RDKit✔️✔️:   if (denom == 0.0) {
    // RDKit✔️✔️:     return 1.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return x / denom;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    range_check(alpha)?;
    range_check(beta)?;
    sparse_same_length(a, b)?;
    tversky_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
        alpha,
        beta,
    )
}

fn range_check(value: f64) -> Result<(), FingerprintError> {
    if !(0.0..=1.0).contains(&value) {
        return Err(FingerprintError::RangeError { value });
    }
    Ok(())
}

fn tversky_core(common: u32, y: u32, z: u32, a: f64, b: f64) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    if y == 0 || z == 0 {
        return Ok(0.0);
    }
    let denom = a * f64::from(y) + b * f64::from(z) + (1.0 - a - b) * x;
    if denom == 0.0 { Ok(1.0) } else { Ok(x / denom) }
}

/// `CosineSimilarity` (BitOps.cpp:320-334).
pub fn cosine(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double CosineSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y * z > 0.0) {
    // RDKit✔️✔️:     return x / sqrt(y * z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    cosine_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `CosineSimilarity` over the sparse bit vectors (same source body).
pub fn sparse_cosine(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double CosineSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y * z > 0.0) {
    // RDKit✔️✔️:     return x / sqrt(y * z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    cosine_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn cosine_core(common: u32, y: u32, z: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    let y = f64::from(y);
    let z = f64::from(z);
    if y * z > 0.0 {
        Ok(x / (y * z).sqrt())
    } else {
        Ok(0.0)
    }
}

/// `KulczynskiSimilarity` (BitOps.cpp:336-350).
pub fn kulczynski(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double KulczynskiSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y * z > 0.0) {
    // RDKit✔️✔️:     return x * (y + z) / (2 * y * z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    kulczynski_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `KulczynskiSimilarity` over the sparse bit vectors (same source body).
pub fn sparse_kulczynski(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double KulczynskiSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y * z > 0.0) {
    // RDKit✔️✔️:     return x * (y + z) / (2 * y * z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    kulczynski_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn kulczynski_core(common: u32, y: u32, z: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    let y = f64::from(y);
    let z = f64::from(z);
    if y * z > 0.0 {
        Ok(x * (y + z) / (2.0 * y * z))
    } else {
        Ok(0.0)
    }
}

/// `DiceSimilarity` (BitOps.cpp:352-366).
pub fn dice(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double DiceSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y + z > 0.0) {
    // RDKit✔️✔️:     return 2 * x / (y + z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    dice_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `DiceSimilarity` over the sparse bit vectors (same source body).
pub fn sparse_dice(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double DiceSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y + z > 0.0) {
    // RDKit✔️✔️:     return 2 * x / (y + z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    dice_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn dice_core(common: u32, y: u32, z: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    let y = f64::from(y);
    let z = f64::from(z);
    if y + z > 0.0 {
        Ok(2.0 * x / (y + z))
    } else {
        Ok(0.0)
    }
}

/// `SokalSimilarity` (BitOps.cpp:368-381). The zero guard tests the
/// integer on-counts (`y == 0 || z == 0`), distinct from the
/// floating-point product guards of Cosine/Kulczynski.
pub fn sokal(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double SokalSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   auto y = bv1.getNumOnBits();
    // RDKit✔️✔️:   auto z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y == 0 || z == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return x / (2. * y + 2. * z - 3. * x);
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    sokal_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `SokalSimilarity` over the sparse bit vectors (same source body).
pub fn sparse_sokal(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double SokalSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   auto y = bv1.getNumOnBits();
    // RDKit✔️✔️:   auto z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y == 0 || z == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return x / (2. * y + 2. * z - 3. * x);
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    sokal_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn sokal_core(common: u32, y: u32, z: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    if y == 0 || z == 0 {
        return Ok(0.0);
    }
    Ok(x / (2.0 * f64::from(y) + 2.0 * f64::from(z) - 3.0 * x))
}

/// `McConnaugheySimilarity` (BitOps.cpp:383-397).
pub fn mcconnaughey(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double McConnaugheySimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y * z > 0.0) {
    // RDKit✔️✔️:     return (x * (y + z) - (y * z)) / (y * z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    mcconnaughey_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `McConnaugheySimilarity` over the sparse bit vectors (same source body).
pub fn sparse_mcconnaughey(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double McConnaugheySimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y * z > 0.0) {
    // RDKit✔️✔️:     return (x * (y + z) - (y * z)) / (y * z);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    mcconnaughey_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn mcconnaughey_core(common: u32, y: u32, z: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    let y = f64::from(y);
    let z = f64::from(z);
    if y * z > 0.0 {
        Ok((x * (y + z) - (y * z)) / (y * z))
    } else {
        Ok(0.0)
    }
}

/// `AsymmetricSimilarity` (BitOps.cpp:399-414): non-commutative —
/// `common / min(y, z)`.
pub fn asymmetric(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double AsymmetricSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   double min = std::min(y, z);
    // RDKit✔️✔️:   if (min > 0.0) {
    // RDKit✔️✔️:     return x / min;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    asymmetric_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `AsymmetricSimilarity` over the sparse bit vectors (same source body).
pub fn sparse_asymmetric(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double AsymmetricSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   double min = std::min(y, z);
    // RDKit✔️✔️:   if (min > 0.0) {
    // RDKit✔️✔️:     return x / min;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    asymmetric_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn asymmetric_core(common: u32, y: u32, z: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    let min = f64::from(y.min(z));
    if min > 0.0 { Ok(x / min) } else { Ok(0.0) }
}

/// `BraunBlanquetSimilarity` (BitOps.cpp:416-431): `common / max(y, z)`.
pub fn braun_blanquet(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double BraunBlanquetSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   double max = std::max(y, z);
    // RDKit✔️✔️:   if (max > 0.0) {
    // RDKit✔️✔️:     return x / max;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    braun_blanquet_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `BraunBlanquetSimilarity` over the sparse bit vectors (same body).
pub fn sparse_braun_blanquet(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double BraunBlanquetSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double y = bv1.getNumOnBits();
    // RDKit✔️✔️:   double z = bv2.getNumOnBits();
    // RDKit✔️✔️:   double max = std::max(y, z);
    // RDKit✔️✔️:   if (max > 0.0) {
    // RDKit✔️✔️:     return x / max;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    braun_blanquet_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn braun_blanquet_core(common: u32, y: u32, z: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    let max = f64::from(y.max(z));
    if max > 0.0 { Ok(x / max) } else { Ok(0.0) }
}

/// `RusselSimilarity` (BitOps.cpp:433-441): `common / bv1.getNumBits()` —
/// the divisor is the vector **length** (not the on-count; oracle:
/// a{1,3} vs b{3,5} over 8 bits = 1/8 = 0.125, and an empty 8-bit pair
/// gives 0/8 = 0.0, not NaN). Non-commutative in principle; the
/// length-equality precondition keeps both sides' divisor equal. A
/// zero-length dense pair raises the `NumOnBitsInCommon` bitmap
/// precondition (oracle-verified); a zero-length sparse pair computes
/// 0/0 = NaN through IEEE double arithmetic.
pub fn russel(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double RusselSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   return x / bv1.getNumBits();
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    let x = f64::from(num_on_bits_in_common(a, b)?);
    Ok(x / f64::from(a.n_bits()))
}

/// `RusselSimilarity` over the sparse bit vectors (same source body; the
/// sparse `NumOnBitsInCommon` template has no bitmap precondition, so a
/// zero-length sparse pair computes 0/0 = NaN directly).
pub fn sparse_russel(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double RusselSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   return x / bv1.getNumBits();
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    let x = f64::from(sparse_num_on_bits_in_common(a, b)?);
    Ok(x / f64::from(a.n_bits()))
}

/// `RogotGoldbergSimilarity` (BitOps.cpp:443-467) with the source's exact
/// branch order (`y == 0 || z == 0 -> 0.0`, then `x == l || d == l ->
/// 1.0`, then zero-denominator -> 0.0).
pub fn rogot_goldberg(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double RogotGoldbergSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   auto y = bv1.getNumOnBits();
    // RDKit✔️✔️:   auto z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y == 0 || z == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double l = bv1.getNumBits();
    // RDKit✔️✔️:   double d = l - y - z + x;
    // RDKit✔️✔️:   double denom1 = y + z;
    // RDKit✔️✔️:   double denom2 = 2 * l - y - z;
    // RDKit✔️✔️:   if ((x == l) || (d == l)) {
    // RDKit✔️✔️:     return 1.0;
    // RDKit✔️✔️:   } else if (denom1 == 0 || denom2 == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return (x / (y + z) + (d) / (2 * l - y - z));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    rogot_goldberg_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
        a.n_bits(),
    )
}

/// `RogotGoldbergSimilarity` over the sparse bit vectors (same body).
pub fn sparse_rogot_goldberg(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double RogotGoldbergSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double x = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   auto y = bv1.getNumOnBits();
    // RDKit✔️✔️:   auto z = bv2.getNumOnBits();
    // RDKit✔️✔️:   if (y == 0 || z == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double l = bv1.getNumBits();
    // RDKit✔️✔️:   double d = l - y - z + x;
    // RDKit✔️✔️:   double denom1 = y + z;
    // RDKit✔️✔️:   double denom2 = 2 * l - y - z;
    // RDKit✔️✔️:   if ((x == l) || (d == l)) {
    // RDKit✔️✔️:     return 1.0;
    // RDKit✔️✔️:   } else if (denom1 == 0 || denom2 == 0) {
    // RDKit✔️✔️:     return 0.0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return (x / (y + z) + (d) / (2 * l - y - z));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    rogot_goldberg_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
        a.n_bits(),
    )
}

fn rogot_goldberg_core(common: u32, y: u32, z: u32, l: u32) -> Result<f64, FingerprintError> {
    let x = f64::from(common);
    if y == 0 || z == 0 {
        return Ok(0.0);
    }
    let y = f64::from(y);
    let z = f64::from(z);
    let l = f64::from(l);
    let d = l - y - z + x;
    let denom1 = y + z;
    let denom2 = 2.0 * l - y - z;
    if (x == l) || (d == l) {
        Ok(1.0)
    } else if denom1 == 0.0 || denom2 == 0.0 {
        Ok(0.0)
    } else {
        Ok(x / (y + z) + d / (2.0 * l - y - z))
    }
}

/// `OnBitSimilarity` (BitOps.cpp:481-495): `common / (a|b).on_count`.
pub fn on_bit(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double OnBitSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double num = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double denom = (bv1 | bv2).getNumOnBits();
    // RDKit✔️✔️:   if (denom > 0) {
    // RDKit✔️✔️:     return num / denom;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    on_bit_core(num_on_bits_in_common(a, b)?, a.or(b)?.num_on_bits())
}

/// `OnBitSimilarity` over the sparse bit vectors (same body; the sparse
/// union retains out-of-length members, which the source counts in the
/// denominator exactly like the raw set union).
pub fn sparse_on_bit(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double OnBitSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double num = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   double denom = (bv1 | bv2).getNumOnBits();
    // RDKit✔️✔️:   if (denom > 0) {
    // RDKit✔️✔️:     return num / denom;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    on_bit_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.union(b).num_on_bits(),
    )
}

fn on_bit_core(common: u32, denom: u32) -> Result<f64, FingerprintError> {
    let num = f64::from(common);
    let denom = f64::from(denom);
    if denom > 0.0 {
        Ok(num / denom)
    } else {
        Ok(0.0)
    }
}

/// `AllBitSimilarity` (BitOps.cpp:538-545): `NumBitsInCommon / n_bits`
/// (Manhattan similarity over on and off bits). A zero-length pair
/// divides 0/0 through IEEE double arithmetic and yields NaN exactly
/// like the pinned build (oracle-verified: `nan`).
pub fn all_bit(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double AllBitSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return double(NumBitsInCommon(bv1, bv2)) / bv1.getNumBits();
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    Ok(f64::from(num_bits_in_common(a, b)?) / f64::from(a.n_bits()))
}

/// `AllBitSimilarity` over the sparse bit vectors (same body; the count
/// difference uses the source's wrapping unsigned arithmetic, so retained
/// out-of-length members are included exactly like the raw sets).
pub fn sparse_all_bit(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: double AllBitSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return double(NumBitsInCommon(bv1, bv2)) / bv1.getNumBits();
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    Ok(f64::from(sparse_num_bits_in_common(a, b)?) / f64::from(a.n_bits()))
}

/// `OnBitProjSimilarity` (BitOps.cpp:620-632): `[common / y, common /
/// z]`, or `[0.0, 0.0]` when the common count is zero. Note the source
/// divides by the on-counts directly with no zero guard beyond `num`.
pub fn on_bit_proj_similarity(
    a: &Fingerprint,
    b: &Fingerprint,
) -> Result<[f64; 2], FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: DoubleVect OnBitProjSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   DoubleVect res(2, 0.0);
    // RDKit✔️✔️:   double num = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   if (num) {
    // RDKit✔️✔️:     res[0] = num / bv1.getNumOnBits();
    // RDKit✔️✔️:     res[1] = num / bv2.getNumOnBits();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    on_bit_proj_core(
        num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

/// `OnBitProjSimilarity` over the sparse bit vectors (same body).
pub fn sparse_on_bit_proj_similarity(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<[f64; 2], FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: DoubleVect OnBitProjSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   DoubleVect res(2, 0.0);
    // RDKit✔️✔️:   double num = NumOnBitsInCommon(bv1, bv2);
    // RDKit✔️✔️:   if (num) {
    // RDKit✔️✔️:     res[0] = num / bv1.getNumOnBits();
    // RDKit✔️✔️:     res[1] = num / bv2.getNumOnBits();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    on_bit_proj_core(
        sparse_num_on_bits_in_common(a, b)?,
        a.num_on_bits(),
        b.num_on_bits(),
    )
}

fn on_bit_proj_core(common: u32, y: u32, z: u32) -> Result<[f64; 2], FingerprintError> {
    let mut res = [0.0f64; 2];
    let num = f64::from(common);
    if num != 0.0 {
        res[0] = num / f64::from(y);
        res[1] = num / f64::from(z);
    }
    Ok(res)
}

/// `OffBitProjSimilarity` (BitOps.cpp:659-671):
/// `[(a|b).off_count / a.off_count, (a|b).off_count / b.off_count]`, or
/// `[0.0, 0.0]` when the union's off-count is zero.
pub fn off_bit_proj_similarity(
    a: &Fingerprint,
    b: &Fingerprint,
) -> Result<[f64; 2], FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: DoubleVect OffBitProjSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   DoubleVect res(2, 0.0);
    // RDKit✔️✔️:   double num = (bv1 | bv2).getNumOffBits();
    // RDKit✔️✔️:   if (num) {
    // RDKit✔️✔️:     res[0] = num / bv1.getNumOffBits();
    // RDKit✔️✔️:     res[1] = num / bv2.getNumOffBits();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    dense_same_length(a, b)?;
    off_bit_proj_core(a.or(b)?.num_off_bits(), a.num_off_bits(), b.num_off_bits())
}

/// `OffBitProjSimilarity` over the sparse bit vectors (same body; the
/// union's off-count uses the declared length, ignoring retained
/// out-of-length members exactly like the raw sets).
pub fn sparse_off_bit_proj_similarity(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
) -> Result<[f64; 2], FingerprintError> {
    // RDKit✔️✔️: template <typename T1, typename T2>
    // RDKit✔️✔️: DoubleVect OffBitProjSimilarity(const T1 &bv1, const T2 &bv2) {
    // RDKit✔️✔️:   if (bv1.getNumBits() != bv2.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("BitVects must be same length");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   DoubleVect res(2, 0.0);
    // RDKit✔️✔️:   double num = (bv1 | bv2).getNumOffBits();
    // RDKit✔️✔️:   if (num) {
    // RDKit✔️✔️:     res[0] = num / bv1.getNumOffBits();
    // RDKit✔️✔️:     res[1] = num / bv2.getNumOffBits();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    sparse_same_length(a, b)?;
    off_bit_proj_core(
        a.union(b).num_off_bits(),
        a.num_off_bits(),
        b.num_off_bits(),
    )
}

fn off_bit_proj_core(num: u32, y_off: u32, z_off: u32) -> Result<[f64; 2], FingerprintError> {
    let mut res = [0.0f64; 2];
    let num = f64::from(num);
    if num != 0.0 {
        res[0] = num / f64::from(y_off);
        res[1] = num / f64::from(z_off);
    }
    Ok(res)
}

/// `SimilarityWrapper` for a two-argument metric (BitOps.h:30-50): the
/// longer operand is folded by the floor length ratio before the metric
/// runs, and `return_distance` computes `1.0 - res`.
///
/// Source-defined edges reproduced: a non-divisible ratio folds by the
/// floor factor (often leaving lengths unequal, so the metric then raises
/// the length error — oracle: 12 vs 8 -> ValueError), and equal lengths
/// call the metric directly (0 vs 0 -> the metric's own behavior).
///
/// Resolved boundary (Supervisor correction queue C4, 2026-09-24):
/// exactly one zero-length operand makes the source compute `numBits / 0`
/// for the fold factor — an actually executed unsigned division by zero
/// (oracle probe: the pinned build dies with SIGFPE, uncatchable as an
/// exception). The Rust port rejects that executed division with
/// `InvalidArguments` carrying the stable fold-division reason: CK safe
/// rejection of undefined source execution, not an RDKit exception and
/// not a full-input parity claim. Zero/zero lengths take the source's
/// equal-length metric branch (no division executes).
pub fn similarity_wrapper(
    a: &Fingerprint,
    b: &Fingerprint,
    metric: fn(&Fingerprint, &Fingerprint) -> Result<f64, FingerprintError>,
    return_distance: bool,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: double SimilarityWrapper(const T &bv1, const T &bv2,
    // RDKit✔️✔️:                          double (*metric)(const T &, const T &),
    // RDKit✔️✔️:                          bool returnDistance = false) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   if (bv1.getNumBits() > bv2.getNumBits()) {
    // RDKit✔️✔️:     T *bv1tmp = FoldFingerprint(bv1, bv1.getNumBits() / bv2.getNumBits());
    // RDKit✔️✔️:     res = metric(*bv1tmp, bv2);
    // RDKit✔️✔️:     delete bv1tmp;
    // RDKit✔️✔️:   } else if (bv2.getNumBits() > bv1.getNumBits()) {
    // RDKit✔️✔️:     T *bv2tmp = FoldFingerprint(bv2, bv2.getNumBits() / bv1.getNumBits());
    // RDKit✔️✔️:     res = metric(bv1, *bv2tmp);
    // RDKit✔️✔️:     delete bv2tmp;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = metric(bv1, bv2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (returnDistance) {
    // RDKit✔️✔️:     res = 1.0 - res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let res = if a.n_bits() > b.n_bits() {
        if b.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint(a, a.n_bits() / b.n_bits())?;
        metric(&folded, b)?
    } else if b.n_bits() > a.n_bits() {
        if a.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint(b, b.n_bits() / a.n_bits())?;
        metric(a, &folded)?
    } else {
        metric(a, b)?
    };
    Ok(if return_distance { 1.0 - res } else { res })
}

fn zero_length_fold_boundary() -> FingerprintError {
    FingerprintError::InvalidArguments {
        reason: "SimilarityWrapper: exactly one zero-length operand makes the \
                 source fold-factor length division a division by zero; CK \
                 safe rejection of undefined source execution",
    }
}

/// `SimilarityWrapper` over the sparse bit vectors (BitOps.h:30-50),
/// folding the longer operand by the floor length ratio without dense
/// materialization; the same zero-length boundary applies.
pub fn sparse_similarity_wrapper(
    a: &SparseBitFingerprint,
    b: &SparseBitFingerprint,
    metric: fn(&SparseBitFingerprint, &SparseBitFingerprint) -> Result<f64, FingerprintError>,
    return_distance: bool,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: double SimilarityWrapper(const T &bv1, const T &bv2,
    // RDKit✔️✔️:                          double (*metric)(const T &, const T &),
    // RDKit✔️✔️:                          bool returnDistance = false) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   if (bv1.getNumBits() > bv2.getNumBits()) {
    // RDKit✔️✔️:     T *bv1tmp = FoldFingerprint(bv1, bv1.getNumBits() / bv2.getNumBits());
    // RDKit✔️✔️:     res = metric(*bv1tmp, bv2);
    // RDKit✔️✔️:     delete bv1tmp;
    // RDKit✔️✔️:   } else if (bv2.getNumBits() > bv1.getNumBits()) {
    // RDKit✔️✔️:     T *bv2tmp = FoldFingerprint(bv2, bv2.getNumBits() / bv1.getNumBits());
    // RDKit✔️✔️:     res = metric(bv1, *bv2tmp);
    // RDKit✔️✔️:     delete bv2tmp;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = metric(bv1, bv2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (returnDistance) {
    // RDKit✔️✔️:     res = 1.0 - res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let res = if a.n_bits() > b.n_bits() {
        if b.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint_sparse(a, a.n_bits() / b.n_bits())?;
        metric(&folded, b)?
    } else if b.n_bits() > a.n_bits() {
        if a.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint_sparse(b, b.n_bits() / a.n_bits())?;
        metric(a, &folded)?
    } else {
        metric(a, b)?
    };
    Ok(if return_distance { 1.0 - res } else { res })
}

/// `SimilarityWrapper` overload for a parameterized metric (BitOps.h:
/// 52-72): identical fold-and-distance structure with `a`/`b` forwarded
/// unchanged; the metric's own parameter validation (Tversky's
/// `RANGE_CHECK`s) runs inside the metric after folding, exactly like the
/// source.
pub fn similarity_wrapper_tversky(
    x: &Fingerprint,
    y: &Fingerprint,
    alpha: f64,
    beta: f64,
    metric: fn(&Fingerprint, &Fingerprint, f64, f64) -> Result<f64, FingerprintError>,
    return_distance: bool,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: double SimilarityWrapper(const T &bv1, const T &bv2, double a, double b,
    // RDKit✔️✔️:                          double (*metric)(const T &, const T &, double, double),
    // RDKit✔️✔️:                          bool returnDistance = false) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   if (bv1.getNumBits() > bv2.getNumBits()) {
    // RDKit✔️✔️:     T *bv1tmp = FoldFingerprint(bv1, bv1.getNumBits() / bv2.getNumBits());
    // RDKit✔️✔️:     res = metric(*bv1tmp, bv2, a, b);
    // RDKit✔️✔️:     delete bv1tmp;
    // RDKit✔️✔️:   } else if (bv2.getNumBits() > bv1.getNumBits()) {
    // RDKit✔️✔️:     T *bv2tmp = FoldFingerprint(bv2, bv2.getNumBits() / bv1.getNumBits());
    // RDKit✔️✔️:     res = metric(bv1, *bv2tmp, a, b);
    // RDKit✔️✔️:     delete bv2tmp;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = metric(bv1, bv2, a, b);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (returnDistance) {
    // RDKit✔️✔️:     res = 1.0 - res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let res = if x.n_bits() > y.n_bits() {
        if y.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint(x, x.n_bits() / y.n_bits())?;
        metric(&folded, y, alpha, beta)?
    } else if y.n_bits() > x.n_bits() {
        if x.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint(y, y.n_bits() / x.n_bits())?;
        metric(x, &folded, alpha, beta)?
    } else {
        metric(x, y, alpha, beta)?
    };
    Ok(if return_distance { 1.0 - res } else { res })
}

/// Sparse twin of the parameterized wrapper (same source body).
pub fn sparse_similarity_wrapper_tversky(
    x: &SparseBitFingerprint,
    y: &SparseBitFingerprint,
    alpha: f64,
    beta: f64,
    metric: fn(
        &SparseBitFingerprint,
        &SparseBitFingerprint,
        f64,
        f64,
    ) -> Result<f64, FingerprintError>,
    return_distance: bool,
) -> Result<f64, FingerprintError> {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: double SimilarityWrapper(const T &bv1, const T &bv2, double a, double b,
    // RDKit✔️✔️:                          double (*metric)(const T &, const T &, double, double),
    // RDKit✔️✔️:                          bool returnDistance = false) {
    // RDKit✔️✔️:   double res = 0.0;
    // RDKit✔️✔️:   if (bv1.getNumBits() > bv2.getNumBits()) {
    // RDKit✔️✔️:     T *bv1tmp = FoldFingerprint(bv1, bv1.getNumBits() / bv2.getNumBits());
    // RDKit✔️✔️:     res = metric(*bv1tmp, bv2, a, b);
    // RDKit✔️✔️:     delete bv1tmp;
    // RDKit✔️✔️:   } else if (bv2.getNumBits() > bv1.getNumBits()) {
    // RDKit✔️✔️:     T *bv2tmp = FoldFingerprint(bv2, bv2.getNumBits() / bv1.getNumBits());
    // RDKit✔️✔️:     res = metric(bv1, *bv2tmp, a, b);
    // RDKit✔️✔️:     delete bv2tmp;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = metric(bv1, bv2, a, b);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (returnDistance) {
    // RDKit✔️✔️:     res = 1.0 - res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    let res = if x.n_bits() > y.n_bits() {
        if y.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint_sparse(x, x.n_bits() / y.n_bits())?;
        metric(&folded, y, alpha, beta)?
    } else if y.n_bits() > x.n_bits() {
        if x.n_bits() == 0 {
            return Err(zero_length_fold_boundary());
        }
        let folded = crate::folding::fold_fingerprint_sparse(y, y.n_bits() / x.n_bits())?;
        metric(x, &folded, alpha, beta)?
    } else {
        metric(x, y, alpha, beta)?
    };
    Ok(if return_distance { 1.0 - res } else { res })
}

macro_rules! calc_vect_params_impl {
    ($fname:ident, $vec:ty) => {
        /// `calcVectParams` (SparseIntVect.h:434-495): shared accumulation for the
        /// count-valued similarity metrics — `v1_sum`/`v2_sum` are the sums of
        /// **absolute** entry values and `and_sum` the sum of per-common-key
        /// minima of absolute values.
        ///
        /// The source's merge walk accumulates into `double`s, so no integer
        /// overflow exists on that path; the only signed-i32 operation is
        /// `abs(entry)`, which is C++-undefined for `i32::MIN` and is reported
        /// per the Supervisor resolution of 2026-09-24 (CK safety boundary).
        pub fn $fname(v1: &$vec, v2: &$vec) -> Result<(f64, f64, f64), FingerprintError> {
            // RDKit✔️✔️: template <typename IndexType>
            // RDKit✔️✔️: void calcVectParams(const SparseIntVect<IndexType> &v1,
            // RDKit✔️✔️:                     const SparseIntVect<IndexType> &v2, double &v1Sum,
            // RDKit✔️✔️:                     double &v2Sum, double &andSum) {
            // RDKit✔️✔️:   if (v1.getLength() != v2.getLength()) {
            // RDKit✔️✔️:     throw ValueErrorException("SparseIntVect size mismatch");
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   v1Sum = v2Sum = andSum = 0.0;
            // RDKit✔️✔️:   // we're doing : (v1&v2).getTotalVal(), but w/o generating
            // RDKit✔️✔️:   // the other vector:
            // RDKit✔️✔️:   typename SparseIntVect<IndexType>::StorageType::const_iterator iter1, iter2;
            // RDKit✔️✔️:   iter1 = v1.getNonzeroElements().begin();
            // RDKit✔️✔️:   if (iter1 != v1.getNonzeroElements().end()) {
            // RDKit✔️✔️:     v1Sum += abs(iter1->second);
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   iter2 = v2.getNonzeroElements().begin();
            // RDKit✔️✔️:   if (iter2 != v2.getNonzeroElements().end()) {
            // RDKit✔️✔️:     v2Sum += abs(iter2->second);
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   while (iter1 != v1.getNonzeroElements().end()) {
            // RDKit✔️✔️:     while (iter2 != v2.getNonzeroElements().end() &&
            // RDKit✔️✔️:            iter2->first < iter1->first) {
            // RDKit✔️✔️:       ++iter2;
            // RDKit✔️✔️:       if (iter2 != v2.getNonzeroElements().end()) {
            // RDKit✔️✔️:         v2Sum += abs(iter2->second);
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:     if (iter2 != v2.getNonzeroElements().end()) {
            // RDKit✔️✔️:       if (iter2->first == iter1->first) {
            // RDKit✔️✔️:         if (abs(iter2->second) < abs(iter1->second)) {
            // RDKit✔️✔️:           andSum += abs(iter2->second);
            // RDKit✔️✔️:         } else {
            // RDKit✔️✔️:           andSum += abs(iter1->second);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         ++iter2;
            // RDKit✔️✔️:         if (iter2 != v2.getNonzeroElements().end()) {
            // RDKit✔️✔️:           v2Sum += abs(iter2->second);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       ++iter1;
            // RDKit✔️✔️:       if (iter1 != v1.getNonzeroElements().end()) {
            // RDKit✔️✔️:         v1Sum += abs(iter1->second);
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       break;
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   if (iter1 != v1.getNonzeroElements().end()) {
            // RDKit✔️✔️:     ++iter1;
            // RDKit✔️✔️:     while (iter1 != v1.getNonzeroElements().end()) {
            // RDKit✔️✔️:       v1Sum += abs(iter1->second);
            // RDKit✔️✔️:       ++iter1;
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   if (iter2 != v2.getNonzeroElements().end()) {
            // RDKit✔️✔️:     ++iter2;
            // RDKit✔️✔️:     while (iter2 != v2.getNonzeroElements().end()) {
            // RDKit✔️✔️:       v2Sum += abs(iter2->second);
            // RDKit✔️✔️:       ++iter2;
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:   }
            // RDKit✔️✔️: }
            //
            // C3 correction: this is the pinned source's monotonic two-iterator
            // merge itself, not an algebraic restatement. `cur1`/`cur2` hold the
            // element the C++ iterators point at (None == end) and `iter1`/
            // `iter2` are positioned just past them, so every `++iter` maps to
            // one `next()` with the same guarded accumulation, match handling and
            // both tail loops. Each accumulator's addition order is the source
            // walk order; `abs` is checked exactly at the executed source abs
            // sites. Complexity: the walk advances each iterator monotonically,
            // so it is O(n + m) comparisons over the two sorted maps with O(1)
            // auxiliary state — no per-key tree lookups, no temporary
            // intersection map, no densification.
            if v1.length() != v2.length() {
                return Err(FingerprintError::BitLengthMismatch {
                    left: u64::from(v1.length()),
                    right: u64::from(v2.length()),
                });
            }
            let checked_abs = |v: i32| -> Result<f64, FingerprintError> {
                v.checked_abs()
                    .map(f64::from)
                    .ok_or(FingerprintError::UndefinedArithmetic {
                        site: "SparseIntVect::calcVectParams abs",
                    })
            };
            let mut v1_sum = 0.0;
            let mut v2_sum = 0.0;
            let mut and_sum = 0.0;
            let mut iter1 = v1.nonzero_elements().iter();
            let mut iter2 = v2.nonzero_elements().iter();
            let mut cur1 = iter1.next();
            if let Some((_, value)) = cur1 {
                v1_sum += checked_abs(*value)?;
            }
            let mut cur2 = iter2.next();
            if let Some((_, value)) = cur2 {
                v2_sum += checked_abs(*value)?;
            }
            while let Some((&key1, v1_value)) = cur1 {
                while let Some((&key2, _)) = cur2 {
                    if key2 < key1 {
                        cur2 = iter2.next();
                        if let Some((_, value)) = cur2 {
                            v2_sum += checked_abs(*value)?;
                        }
                    } else {
                        break;
                    }
                }
                let Some((&key2, &v2_value)) = cur2 else {
                    break;
                };
                if key2 == key1 {
                    let a1 = checked_abs(*v1_value)?;
                    let a2 = checked_abs(v2_value)?;
                    and_sum += if a2 < a1 { a2 } else { a1 };
                    cur2 = iter2.next();
                    if let Some((_, value)) = cur2 {
                        v2_sum += checked_abs(*value)?;
                    }
                }
                cur1 = iter1.next();
                if let Some((_, value)) = cur1 {
                    v1_sum += checked_abs(*value)?;
                }
            }
            if cur1.is_some() {
                // source tail: the pre-advance skips the current element, whose
                // value was already accumulated when the walk advanced to it
                cur1 = iter1.next();
                while let Some((_, value)) = cur1 {
                    v1_sum += checked_abs(*value)?;
                    cur1 = iter1.next();
                }
            }
            if cur2.is_some() {
                cur2 = iter2.next();
                while let Some((_, value)) = cur2 {
                    v2_sum += checked_abs(*value)?;
                    cur2 = iter2.next();
                }
            }
            Ok((v1_sum, v2_sum, and_sum))
        }
    };
}

calc_vect_params_impl!(calc_vect_params, crate::SparseCountFingerprint);
calc_vect_params_impl!(calc_vect_params_u32, crate::SparseCountFingerprint32);

macro_rules! sparse_count_dice_impl {
    ($fname:ident, $vec:ty, $calc:ident) => {
        /// Count-valued `DiceSimilarity` (SparseIntVect.h:498-539), including the
        /// bounds early-exit (only when `!return_distance && bounds > 0.0`) and
        /// `returnDistance`, with the source's exact `fabs(denom) < 1e-6`
        /// thresholds.
        pub fn $fname(
            v1: &$vec,
            v2: &$vec,
            return_distance: bool,
            bounds: f64,
        ) -> Result<f64, FingerprintError> {
            // RDKit✔️✔️: template <typename IndexType>
            // RDKit✔️✔️: double DiceSimilarity(const SparseIntVect<IndexType> &v1,
            // RDKit✔️✔️:                       const SparseIntVect<IndexType> &v2,
            // RDKit✔️✔️:                       bool returnDistance = false, double bounds = 0.0) {
            // RDKit✔️✔️:   if (v1.getLength() != v2.getLength()) {
            // RDKit✔️✔️:     throw ValueErrorException("SparseIntVect size mismatch");
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   double v1Sum = 0.0;
            // RDKit✔️✔️:   double v2Sum = 0.0;
            // RDKit✔️✔️:   if (!returnDistance && bounds > 0.0) {
            // RDKit✔️✔️:     v1Sum = v1.getTotalVal(true);
            // RDKit✔️✔️:     v2Sum = v2.getTotalVal(true);
            // RDKit✔️✔️:     double denom = v1Sum + v2Sum;
            // RDKit✔️✔️:     if (fabs(denom) < 1e-6) {
            // RDKit✔️✔️:       // no need to worry about returnDistance here
            // RDKit✔️✔️:       return 0.0;
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:     double minV = v1Sum < v2Sum ? v1Sum : v2Sum;
            // RDKit✔️✔️:     if (2. * minV / denom < bounds) {
            // RDKit✔️✔️:       return 0.0;
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:     v1Sum = 0.0;
            // RDKit✔️✔️:     v2Sum = 0.0;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   double numer = 0.0;
            // RDKit✔️✔️:   calcVectParams(v1, v2, v1Sum, v2Sum, numer);
            // RDKit✔️✔️:   double denom = v1Sum + v2Sum;
            // RDKit✔️✔️:   double sim;
            // RDKit✔️✔️:   if (fabs(denom) < 1e-6) {
            // RDKit✔️✔️:     sim = 0.0;
            // RDKit✔️✔️:   } else {
            // RDKit✔️✔️:     sim = 2. * numer / denom;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   if (returnDistance) {
            // RDKit✔️✔️:     sim = 1. - sim;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   return sim;
            // RDKit✔️✔️: }
            if v1.length() != v2.length() {
                return Err(FingerprintError::BitLengthMismatch {
                    left: u64::from(v1.length()),
                    right: u64::from(v2.length()),
                });
            }
            if !return_distance && bounds > 0.0 {
                let v1_sum = f64::from(v1.total_value(true)?);
                let v2_sum = f64::from(v2.total_value(true)?);
                let denom = v1_sum + v2_sum;
                if denom.abs() < 1e-6 {
                    return Ok(0.0);
                }
                let min_v = v1_sum.min(v2_sum);
                if 2.0 * min_v / denom < bounds {
                    return Ok(0.0);
                }
            }
            let (v1_sum, v2_sum, numer) = $calc(v1, v2)?;
            let denom = v1_sum + v2_sum;
            let sim = if denom.abs() < 1e-6 {
                0.0
            } else {
                2.0 * numer / denom
            };
            Ok(if return_distance { 1.0 - sim } else { sim })
        }
    };
}

sparse_count_dice_impl!(
    sparse_count_dice,
    crate::SparseCountFingerprint,
    calc_vect_params
);
sparse_count_dice_impl!(
    sparse_count_dice_u32,
    crate::SparseCountFingerprint32,
    calc_vect_params_u32
);

macro_rules! sparse_count_tversky_impl {
    ($fname:ident, $vec:ty, $calc:ident) => {
        /// Count-valued `TverskySimilarity` (SparseIntVect.h:541-568). Unlike the
        /// bit-vector variant, the source performs **no** parameter range checks
        /// here, and `bounds` is explicitly unused (`RDUNUSED_PARAM`); the
        /// denominator threshold is `fabs(denom) < 1e-6`.
        pub fn $fname(
            v1: &$vec,
            v2: &$vec,
            a: f64,
            b: f64,
            return_distance: bool,
            bounds: f64,
        ) -> Result<f64, FingerprintError> {
            // RDKit✔️✔️: template <typename IndexType>
            // RDKit✔️✔️: double TverskySimilarity(const SparseIntVect<IndexType> &v1,
            // RDKit✔️✔️:                          const SparseIntVect<IndexType> &v2, double a, double b,
            // RDKit✔️✔️:                          bool returnDistance = false, double bounds = 0.0) {
            // RDKit✔️✔️:   RDUNUSED_PARAM(bounds);
            // RDKit✔️✔️:   if (v1.getLength() != v2.getLength()) {
            // RDKit✔️✔️:     throw ValueErrorException("SparseIntVect size mismatch");
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   double v1Sum = 0.0;
            // RDKit✔️✔️:   double v2Sum = 0.0;
            // RDKit✔️✔️:   double andSum = 0.0;
            // RDKit✔️✔️:   calcVectParams(v1, v2, v1Sum, v2Sum, andSum);
            // RDKit✔️✔️:   double denom = a * v1Sum + b * v2Sum + (1 - a - b) * andSum;
            // RDKit✔️✔️:   double sim;
            // RDKit✔️✔️:   if (fabs(denom) < 1e-6) {
            // RDKit✔️✔️:     sim = 0.0;
            // RDKit✔️✔️:   } else {
            // RDKit✔️✔️:     sim = andSum / denom;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   if (returnDistance) {
            // RDKit✔️✔️:     sim = 1. - sim;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   return sim;
            // RDKit✔️✔️: }
            let _ = bounds;
            if v1.length() != v2.length() {
                return Err(FingerprintError::BitLengthMismatch {
                    left: u64::from(v1.length()),
                    right: u64::from(v2.length()),
                });
            }
            let (v1_sum, v2_sum, and_sum) = $calc(v1, v2)?;
            let denom = a * v1_sum + b * v2_sum + (1.0 - a - b) * and_sum;
            let sim = if denom.abs() < 1e-6 {
                0.0
            } else {
                and_sum / denom
            };
            Ok(if return_distance { 1.0 - sim } else { sim })
        }
    };
}

sparse_count_tversky_impl!(
    sparse_count_tversky,
    crate::SparseCountFingerprint,
    calc_vect_params
);
sparse_count_tversky_impl!(
    sparse_count_tversky_u32,
    crate::SparseCountFingerprint32,
    calc_vect_params_u32
);

macro_rules! sparse_count_tanimoto_impl {
    ($fname:ident, $vec:ty, $tversky:ident) => {
        /// Count-valued `TanimotoSimilarity` (SparseIntVect.h:570-575): defined
        /// by the source as exactly `TverskySimilarity(v1, v2, 1.0, 1.0, ...)`
        /// with the flags forwarded.
        pub fn $fname(
            v1: &$vec,
            v2: &$vec,
            return_distance: bool,
            bounds: f64,
        ) -> Result<f64, FingerprintError> {
            // RDKit✔️✔️: template <typename IndexType>
            // RDKit✔️✔️: double TanimotoSimilarity(const SparseIntVect<IndexType> &v1,
            // RDKit✔️✔️:                           const SparseIntVect<IndexType> &v2,
            // RDKit✔️✔️:                           bool returnDistance = false, double bounds = 0.0) {
            // RDKit✔️✔️:   return TverskySimilarity(v1, v2, 1.0, 1.0, returnDistance, bounds);
            // RDKit✔️✔️: }
            $tversky(v1, v2, 1.0, 1.0, return_distance, bounds)
        }
    };
}

sparse_count_tanimoto_impl!(
    sparse_count_tanimoto,
    crate::SparseCountFingerprint,
    sparse_count_tversky
);
sparse_count_tanimoto_impl!(
    sparse_count_tanimoto_u32,
    crate::SparseCountFingerprint32,
    sparse_count_tversky_u32
);

/// C2 correction: the source's unsigned Tanimoto total
/// (`bv1.getNumOnBits() + bv2.getNumOnBits()`, BitOps.cpp:291) wraps
/// modulo 2^32; explicit wrapping keeps the behavior independent of the
/// build profile's overflow checks. Private helper actually called by
/// both production Tanimoto overloads.
const fn tanimoto_totals(y: u32, z: u32) -> u32 {
    y.wrapping_add(z)
}

/// C2 correction: the source's unsigned denominator
/// (`total - common`, BitOps.cpp:296). See [`tanimoto_totals`].
const fn tanimoto_denominator(total: u32, common: u32) -> u32 {
    total.wrapping_sub(common)
}

#[cfg(test)]
mod c2_boundary_tests {
    //! Bounded owning-module tests for the Tanimoto unsigned-u32 sites
    //! (C2), exercised directly through the production helpers.
    use super::{tanimoto_denominator, tanimoto_totals};

    #[test]
    fn tanimoto_total_wraps_at_u32_max() {
        assert_eq!(tanimoto_totals(u32::MAX, 1), 0);
        assert_eq!(tanimoto_totals(u32::MAX, 0), u32::MAX);
        assert_eq!(tanimoto_totals(1 << 31, 1 << 31), 0);
        assert_eq!(tanimoto_totals(3, 3), 6);
    }

    #[test]
    fn tanimoto_denominator_wraps_when_common_exceeds_wrapped_total() {
        assert_eq!(tanimoto_denominator(0, 1), u32::MAX);
        assert_eq!(tanimoto_denominator(6, 2), 4);
        assert_eq!(tanimoto_denominator(u32::MAX, u32::MAX), 0);
    }
}
