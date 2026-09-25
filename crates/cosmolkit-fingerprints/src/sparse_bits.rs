//! Sparse fingerprint bit vector: source-backed port of RDKit
//! `Code/DataStructs/SparseBitVect.{h,cpp}`.
//!
//! Carrier fidelity (F02): the source stores on bits in a `std::set<int>`
//! (`IntSet` typedef, SparseBitVect.h:22) with `unsigned int d_size`.
//! `setBit(unsigned)`/`getBit(unsigned)` validate the unsigned index with
//! `checkIndex` and then insert/look up the value converted to `int`: on
//! the pinned ABI (gcc x86-64, oracle-verified against RDKit 2026.03.1)
//! indices >= 2^31 wrap to negative ints (2^31+1 -> -2147483647,
//! u32::MAX -> -1). `getOnBits` returns the raw `int` values in the set's
//! signed ascending order, so high indices sort first and appear negative.
//! The Rust port reproduces this with a `BTreeSet<i32>` carrier and the
//! same wrapping `u32 -> i32` key conversion at every source boundary.

use std::collections::BTreeSet;

use crate::FingerprintError;

/// Sparse bit vector (`SparseBitVect` equivalent).
#[derive(Debug, Clone, Eq)]
pub struct SparseBitFingerprint {
    n_bits: u32,
    on: BTreeSet<i32>,
}

impl PartialEq for SparseBitFingerprint {
    /// Source `operator==` compares only the on-bit sets (`*dp_bits ==
    /// *o.dp_bits`), ignoring `d_size` (SparseBitVect.h:90-92). Pinned by
    /// oracle probe: 4-bit `{1}` == 8-bit `{1}`.
    fn eq(&self, other: &Self) -> bool {
        self.on == other.on
    }
}

/// Source `unsigned int -> int` key conversion for `set<int>` storage:
/// wrapping reinterpret, verified against the pinned build's observable
/// behavior (oracle probes in tests).
const fn to_key(bit: u32) -> i32 {
    bit as i32
}

impl SparseBitFingerprint {
    /// Size constructor.
    ///
    /// Source: `explicit SparseBitVect(unsigned int size)` +
    /// `_initForSize` (SparseBitVect.cpp:39-41, 311-315).
    #[must_use]
    pub fn new(n_bits: u32) -> Self {
        // RDKit✔️✔️: explicit SparseBitVect(unsigned int size) : dp_bits(nullptr), d_size(0) {
        // RDKit✔️✔️:   _initForSize(size);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: void SparseBitVect::_initForSize(unsigned int size) {
        // RDKit✔️✔️:   d_size = size;
        // RDKit✔️✔️:   delete dp_bits;
        // RDKit✔️✔️:   dp_bits = new IntSet;
        // RDKit✔️✔️: };
        Self {
            n_bits,
            on: BTreeSet::new(),
        }
    }

    /// Number of bits.
    ///
    /// Source: `unsigned int getNumBits() const override { return d_size; }`
    /// (SparseBitVect.h:68).
    #[must_use]
    pub const fn n_bits(&self) -> u32 {
        // RDKit✔️✔️: unsigned int getNumBits() const override { return d_size; }
        self.n_bits
    }

    /// Source index check.
    ///
    /// Source: `bool checkIndex(const unsigned int idx) const { return idx <
    /// d_size || (idx == d_size && d_size ==
    /// std::numeric_limits<unsigned int>::max()); }`
    /// (SparseBitVect.h:100-103).
    const fn check_index(&self, idx: u32) -> bool {
        // RDKit✔️✔️: bool checkIndex(const unsigned int idx) const {
        // RDKit✔️✔️:   return idx < d_size || (idx == d_size &&
        // RDKit✔️✔️:                             d_size == std::numeric_limits<unsigned int>::max());
        // RDKit✔️✔️: }
        idx < self.n_bits || (idx == self.n_bits && self.n_bits == u32::MAX)
    }

    fn checked(&self, idx: u32) -> Result<(), FingerprintError> {
        if !self.check_index(idx) {
            return Err(FingerprintError::SparseIndexOutOfRange {
                index: u64::from(idx),
                size: u64::from(self.n_bits),
            });
        }
        Ok(())
    }

    /// Query one bit by unsigned index.
    ///
    /// Source: `bool SparseBitVect::getBit(const unsigned int which)
    /// const` (SparseBitVect.cpp:150-155). The lookup converts `which` to
    /// the signed stored key, wrapping for indices >= 2^31.
    pub fn get_bit(&self, bit: u32) -> Result<bool, FingerprintError> {
        // RDKit✔️✔️: bool SparseBitVect::getBit(const unsigned int which) const {
        // RDKit✔️✔️:   if (!checkIndex(which)) {
        // RDKit✔️✔️:     throw IndexErrorException(which);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return dp_bits->count(which) > 0u;
        // RDKit✔️✔️: }
        self.checked(bit)?;
        Ok(self.on.contains(&to_key(bit)))
    }

    /// Set one bit by unsigned index; returns its original state.
    ///
    /// Source: `bool SparseBitVect::setBit(const unsigned int which)`
    /// (SparseBitVect.cpp:190-199). `dp_bits->insert(which)` converts the
    /// unsigned index to the signed stored key.
    pub fn set_bit(&mut self, bit: u32) -> Result<bool, FingerprintError> {
        // RDKit✔️✔️: bool SparseBitVect::setBit(const unsigned int which) {
        // RDKit✔️✔️:   if (!dp_bits) {
        // RDKit✔️✔️:     throw ValueErrorException("BitVect not properly initialized.");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!checkIndex(which)) {
        // RDKit✔️✔️:     throw IndexErrorException(which);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   auto res = dp_bits->insert(which);
        // RDKit✔️✔️:   return !(res.second);
        // RDKit✔️✔️: }
        //
        // The `!dp_bits` ValueErrorException branch is unreachable in the
        // Rust ownership model: `on` always exists. Recorded, not modeled.
        self.checked(bit)?;
        Ok(!self.on.insert(to_key(bit)))
    }

    /// Unset one bit by unsigned index; returns its original state.
    ///
    /// Source: `bool SparseBitVect::unsetBit(const unsigned int which)`
    /// (SparseBitVect.cpp:227-241).
    pub fn unset_bit(&mut self, bit: u32) -> Result<bool, FingerprintError> {
        // RDKit✔️✔️: bool SparseBitVect::unsetBit(const unsigned int which) {
        // RDKit✔️✔️:   if (!dp_bits) {
        // RDKit✔️✔️:     throw ValueErrorException("BitVect not properly initialized.");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!checkIndex(which)) {
        // RDKit✔️✔️:     throw IndexErrorException(which);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (dp_bits->count(which)) {
        // RDKit✔️✔️:     dp_bits->erase(dp_bits->find(which));
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        self.checked(bit)?;
        Ok(self.on.remove(&to_key(bit)))
    }

    /// Number of on bits.
    ///
    /// Source: `unsigned int getNumOnBits() const override { return
    /// static_cast<unsigned int>(dp_bits->size()); }` (SparseBitVect.h:76-78).
    #[must_use]
    pub fn num_on_bits(&self) -> u32 {
        // RDKit✔️✔️: unsigned int getNumOnBits() const override {
        // RDKit✔️✔️:   return static_cast<unsigned int>(dp_bits->size());
        // RDKit✔️✔️: }
        self.on.len() as u32
    }

    /// Number of off bits.
    ///
    /// Source: `unsigned int getNumOffBits() const override { return
    /// d_size - static_cast<unsigned int>(dp_bits->size()); }`
    /// (SparseBitVect.h:79-81).
    #[must_use]
    pub fn num_off_bits(&self) -> u32 {
        // RDKit✔️✔️: unsigned int getNumOffBits() const override {
        // RDKit✔️✔️:   return d_size - static_cast<unsigned int>(dp_bits->size());
        // RDKit✔️✔️: }
        self.n_bits - self.num_on_bits()
    }

    /// On bits as raw signed storage values, in the set's signed ascending
    /// order (source `IntVect` is `std::vector<int>`).
    ///
    /// Source: `void SparseBitVect::getOnBits(IntVect &v) const`
    /// (SparseBitVect.cpp:252-263): the `set<int>` iterates in signed
    /// order, so indices >= 2^31 appear first as wrapped negative values
    /// (oracle: {5, 2^31+1} -> [-2147483647, 5]; {7, u32::MAX} -> [-1, 7]).
    #[must_use]
    pub fn on_bits(&self) -> Vec<i32> {
        // RDKit✔️✔️: void SparseBitVect::getOnBits(IntVect &v) const {
        // RDKit✔️✔️:   if (!dp_bits) {
        // RDKit✔️✔️:     throw ValueErrorException("BitVect not properly initialized.");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   unsigned int nOn = getNumOnBits();
        // RDKit✔️✔️:   if (!v.empty()) {
        // RDKit✔️✔️:     IntVect().swap(v);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   v.reserve(nOn);
        // RDKit✔️✔️:   v.resize(nOn);
        // RDKit✔️✔️:   std::copy(dp_bits->begin(), dp_bits->end(), v.begin());
        // RDKit✔️✔️: };
        self.on.iter().copied().collect()
    }

    /// Access the on-bit set (source `getBitSet()` raw-storage access).
    #[must_use]
    pub const fn bit_set(&self) -> &BTreeSet<i32> {
        // RDKit✔️✔️: const IntSet *getBitSet() const { return dp_bits; }
        &self.on
    }

    /// Clear all bits.
    ///
    /// Source: `void clearBits() override { dp_bits->clear(); }`
    /// (SparseBitVect.h:86).
    pub fn clear_bits(&mut self) {
        // RDKit✔️✔️: void clearBits() override { dp_bits->clear(); }
        self.on.clear()
    }

    /// Set union.
    ///
    /// Source: `SparseBitVect SparseBitVect::operator|(const SparseBitVect
    /// &other) const` (SparseBitVect.cpp:91-97). The source performs **no**
    /// length check: the result carries the left operand's `d_size` and the
    /// raw set union, so members of `other` beyond the result's own length
    /// are retained (oracle probe: 8|16 with right bit 12 -> result 8 bits,
    /// on [1,2,12]; `GetBit(12)` on that result still raises IndexError).
    /// The earlier fail-closed length rejection was an unapproved
    /// divergence and is removed.
    pub fn union(&self, other: &Self) -> Self {
        // RDKit✔️✔️: SparseBitVect SparseBitVect::operator|(const SparseBitVect &other) const {
        // RDKit✔️✔️:   SparseBitVect ans(d_size);
        // RDKit✔️✔️:   std::set_union(dp_bits->begin(), dp_bits->end(), other.dp_bits->begin(),
        // RDKit✔️✔️:                  other.dp_bits->end(),
        // RDKit✔️✔️:                  std::inserter(*(ans.dp_bits), ans.dp_bits->end()));
        // RDKit✔️✔️:   return ans;
        // RDKit✔️✔️: }
        let mut out = Self::new(self.n_bits);
        out.on.extend(self.on.union(&other.on).copied());
        out
    }

    /// Set intersection.
    ///
    /// Source: `SparseBitVect SparseBitVect::operator&(...) const`
    /// (SparseBitVect.cpp:105-111). No length check; result length is the
    /// left operand's.
    pub fn intersection(&self, other: &Self) -> Self {
        // RDKit✔️✔️: SparseBitVect SparseBitVect::operator&(const SparseBitVect &other) const {
        // RDKit✔️✔️:   SparseBitVect ans(d_size);
        // RDKit✔️✔️:   std::set_intersection(dp_bits->begin(), dp_bits->end(),
        // RDKit✔️✔️:                         other.dp_bits->begin(), other.dp_bits->end(),
        // RDKit✔️✔️:                         std::inserter(*(ans.dp_bits), ans.dp_bits->end()));
        // RDKit✔️✔️:   return ans;
        // RDKit✔️✔️: }
        let mut out = Self::new(self.n_bits);
        out.on.extend(self.on.intersection(&other.on).copied());
        out
    }

    /// Set symmetric difference.
    ///
    /// Source: `SparseBitVect SparseBitVect::operator^(...) const`
    /// (SparseBitVect.cpp:119-125). No length check; result length is the
    /// left operand's; out-of-length members are retained.
    pub fn symmetric_difference(&self, other: &Self) -> Self {
        // RDKit✔️✔️: SparseBitVect SparseBitVect::operator^(const SparseBitVect &other) const {
        // RDKit✔️✔️:   SparseBitVect ans(d_size);
        // RDKit✔️✔️:   std::set_symmetric_difference(
        // RDKit✔️✔️:       dp_bits->begin(), dp_bits->end(), other.dp_bits->begin(),
        // RDKit✔️✔️:       other.dp_bits->end(), std::inserter(*(ans.dp_bits), ans.dp_bits->end()));
        // RDKit✔️✔️:   return (ans);
        // RDKit✔️✔️: }
        let mut out = Self::new(self.n_bits);
        out.on
            .extend(self.on.symmetric_difference(&other.on).copied());
        out
    }

    /// Complement (materializes all off bits below `d_size`).
    ///
    /// Source: `SparseBitVect SparseBitVect::operator~() const`
    /// (SparseBitVect.cpp:133-142). The source loops `i < d_size`, so the
    /// `idx == d_size` allowance at `d_size == u32::MAX` is never visited
    /// by the complement; stored out-of-length members do not affect the
    /// loop because `getBit` is only queried below `d_size`.
    #[must_use]
    pub fn negate(&self) -> Self {
        // RDKit✔️✔️: SparseBitVect SparseBitVect::operator~() const {
        // RDKit✔️✔️:   SparseBitVect ans(d_size);
        // RDKit✔️✔️:   for (unsigned int i = 0; i < d_size; i++) {
        // RDKit✔️✔️:     if (!getBit(i)) {
        // RDKit✔️✔️:       ans.setBit(i);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return (ans);
        // RDKit✔️✔️: }
        let mut out = Self::new(self.n_bits);
        for idx in 0..self.n_bits {
            if !self.on.contains(&to_key(idx)) {
                out.on.insert(to_key(idx));
            }
        }
        out
    }
}
