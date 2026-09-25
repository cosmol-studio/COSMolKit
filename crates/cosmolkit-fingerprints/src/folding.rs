//! Fingerprint folding and value conversions: source-backed port of
//! `BitOps.cpp::FoldFingerprint` over the dense and sparse bit vectors.

use crate::{Fingerprint, FingerprintError, SparseBitFingerprint};

/// `FoldFingerprint<ExplicitBitVect>` (BitOps.cpp:673-690).
///
/// Factor contract (source): `factor == 0 || factor >= n_bits` raises
/// `ValueErrorException("invalid fold factor")`; the result length is
/// `init_size / factor` (floor); every on bit maps to `on_bit % res_size`.
pub fn fold_fingerprint(fp: &Fingerprint, factor: u32) -> Result<Fingerprint, FingerprintError> {
    // RDKit✔️✔️: template <typename T1>
    // RDKit✔️✔️: T1 *FoldFingerprint(const T1 &bv1, unsigned int factor) {
    // RDKit✔️✔️:   if (factor <= 0 || factor >= bv1.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("invalid fold factor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int initSize = bv1.getNumBits();
    // RDKit✔️✔️:   int resSize = initSize / factor;
    // RDKit✔️✔️:   auto *res = new T1(resSize);
    // RDKit✔️✔️:   IntVect onBits;
    // RDKit✔️✔️:   bv1.getOnBits(onBits);
    // RDKit✔️✔️:   for (int &onBit : onBits) {
    // RDKit✔️✔️:     int pos = onBit % resSize;
    // RDKit✔️✔️:     res->setBit(pos);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if factor == 0 || factor >= fp.n_bits() {
        return Err(FingerprintError::InvalidFoldFactor {
            factor,
            n_bits: fp.n_bits(),
        });
    }
    let res_size = fp.n_bits() / factor;
    let mut res = Fingerprint::new(res_size);
    for on_bit in fp.on_bits() {
        res.set_bit(on_bit % res_size)?;
    }
    Ok(res)
}

/// `FoldFingerprint<SparseBitVect>` (BitOps.cpp:673-690) over the sparse
/// bit vector, without dense materialization.
///
/// Signed index arithmetic (pinned ABI): `getOnBits` yields the raw
/// `set<int>` keys, so indices >= 2^31 arrive as negative ints. The source
/// computes `int pos = onBit % resSize` with C++ truncation division
/// (negative dividend -> negative remainder) and calls `setBit(pos)`,
/// whose `unsigned int` conversion wraps; a negative `pos` therefore fails
/// `checkIndex` and raises `IndexErrorException` (oracle probe: folding
/// `{2^31+1, 3}` over size 2^31+4 by factor 2 raises IndexError with the
/// negative remainder -1073741821). The Rust error records the wrapped
/// u32 index rather than the C++ `int`-narrowed message value (same
/// boundary, structured fields).
pub fn fold_fingerprint_sparse(
    fp: &SparseBitFingerprint,
    factor: u32,
) -> Result<SparseBitFingerprint, FingerprintError> {
    // RDKit✔️✔️: template <typename T1>
    // RDKit✔️✔️: T1 *FoldFingerprint(const T1 &bv1, unsigned int factor) {
    // RDKit✔️✔️:   if (factor <= 0 || factor >= bv1.getNumBits()) {
    // RDKit✔️✔️:     throw ValueErrorException("invalid fold factor");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   int initSize = bv1.getNumBits();
    // RDKit✔️✔️:   int resSize = initSize / factor;
    // RDKit✔️✔️:   auto *res = new T1(resSize);
    // RDKit✔️✔️:   IntVect onBits;
    // RDKit✔️✔️:   bv1.getOnBits(onBits);
    // RDKit✔️✔️:   for (int &onBit : onBits) {
    // RDKit✔️✔️:     int pos = onBit % resSize;
    // RDKit✔️✔️:     res->setBit(pos);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    if factor == 0 || factor >= fp.n_bits() {
        return Err(FingerprintError::InvalidFoldFactor {
            factor,
            n_bits: fp.n_bits(),
        });
    }
    // `int resSize = initSize / factor` — the C++ unsigned->int conversion
    // wraps on gcc for values above i32::MAX; the constructor's int->
    // unsigned conversion wraps back, so the stored length equals the u32
    // quotient while the modulo below uses the wrapped i32 exactly like
    // the source.
    let res_size = fp.n_bits() / factor;
    let res_size_i32 = res_size as i32;
    let mut res = SparseBitFingerprint::new(res_size);
    for key in fp.on_bits() {
        let pos = key % res_size_i32;
        res.set_bit(pos as u32)?;
    }
    Ok(res)
}

/// Dense -> sparse conversion: length-preserving reconstruction of the
/// on-bit set through `getOnBits` + `setBit` (the way the source builds a
/// `SparseBitVect` from an `ExplicitBitVect` at call sites; DataStructs
/// defines no single named conversion entry). Dense indices are `u32` and
/// strictly below `n_bits`, so the mapping is lossless by construction.
#[must_use]
pub fn to_sparse_bit(fp: &Fingerprint) -> SparseBitFingerprint {
    let mut out = SparseBitFingerprint::new(fp.n_bits());
    for bit in fp.on_bits() {
        out.set_bit(bit)
            .expect("dense on bits are below n_bits, so sparse set_bit is in range");
    }
    out
}

/// Sparse -> dense conversion: length-preserving reconstruction. The
/// sparse carrier stores signed keys (F02): each key converts back to its
/// `u32` index, and the dense `set_bit` check rejects any member outside
/// the vector length (retained out-of-length members from the source's
/// unchecked set operators) instead of truncating silently.
///
/// `from_lsb_bytes` note (F28 caller evidence): the legacy byte-order
/// constructor's only caller in the legacy tree is the Avalon generator
/// (`legacy-core/properties/avalon_fingerprint.rs:157`), which belongs to
/// the excluded FP-generator family; it is therefore outside this lane's
/// frozen closure on caller evidence, not on the earlier unsupported
/// "depiction ownership" claim (corrected here).
pub fn to_dense(sparse: &SparseBitFingerprint) -> Result<Fingerprint, FingerprintError> {
    let mut out = Fingerprint::new(sparse.n_bits());
    for &key in sparse.bit_set() {
        out.set_bit(key as u32)?;
    }
    Ok(out)
}
