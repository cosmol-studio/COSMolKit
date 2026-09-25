//! Dense fingerprint bit vector: source-backed port of RDKit
//! `Code/DataStructs/ExplicitBitVect.{h,cpp}` over `BitVect.h`.
//!
//! Storage mirrors `boost::dynamic_bitset<>` semantics: LSB-first bits in
//! `u64` words, with the invariant that bits at index >= `n_bits` within the
//! final word are always zero.

use crate::FingerprintError;

/// C2 correction: the source's unsigned cached-count increment
/// (`++d_numOnBits`, ExplicitBitVect.cpp:92) wraps modulo 2^32; explicit
/// wrapping keeps the behavior independent of the build profile's
/// overflow checks. Private helper actually called by production.
const fn cached_count_increment(current: u32) -> u32 {
    current.wrapping_add(1)
}

/// C2 correction: the source's unsigned cached-count decrement
/// (`--d_numOnBits`, ExplicitBitVect.cpp:102). See
/// [`cached_count_increment`].
const fn cached_count_decrement(current: u32) -> u32 {
    current.wrapping_sub(1)
}

/// C2 correction: the source's unsigned off-count subtraction
/// (`d_size - d_numOnBits`, ExplicitBitVect.cpp:182). See
/// [`cached_count_increment`].
const fn off_count(d_size: u32, d_num_on_bits: u32) -> u32 {
    d_size.wrapping_sub(d_num_on_bits)
}

/// Dense bit vector (`ExplicitBitVect` equivalent).
#[derive(Debug, Clone, Eq)]
pub struct Fingerprint {
    words: Vec<u64>,
    n_bits: u32,
    num_on: u32,
}

impl PartialEq for Fingerprint {
    /// Source equality (C1 correction): `ExplicitBitVect::operator==`
    /// compares only the bitsets (`*dp_bits == *o.dp_bits`,
    /// ExplicitBitVect.h:86-91); `dynamic_bitset::==` compares size and
    /// all bits. The cached `d_numOnBits` — which the source deliberately
    /// leaves stale across `clearBits` — takes no part in equality, so
    /// the Rust port compares `n_bits` and `words` only and never
    /// recomputes or clears the cache to make equality pass.
    fn eq(&self, other: &Self) -> bool {
        // RDKit✔️✔️: bool operator==(const ExplicitBitVect &o) const {
        // RDKit✔️✔️:   return *dp_bits == *o.dp_bits;
        // RDKit✔️✔️: }
        self.n_bits == other.n_bits && self.words == other.words
    }
}

impl Fingerprint {
    /// Zero-length vector (source default constructor `ExplicitBitVect()`).
    #[must_use]
    pub const fn empty() -> Self {
        Self {
            words: Vec::new(),
            n_bits: 0,
            num_on: 0,
        }
    }

    /// Size-only constructor.
    ///
    /// Source: `ExplicitBitVect::ExplicitBitVect(unsigned int size)` +
    /// `_initForSize` (ExplicitBitVect.cpp:33-36, 199-203).
    #[must_use]
    pub fn new(n_bits: u32) -> Self {
        // RDKit✔️✔️: explicit ExplicitBitVect(unsigned int size)
        // RDKit✔️✔️:     : dp_bits(nullptr), d_size(0), d_numOnBits(0) {
        // RDKit✔️✔️:   _initForSize(size);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: void ExplicitBitVect::_initForSize(unsigned int size) {
        // RDKit✔️✔️:   d_size = size;
        // RDKit✔️✔️:   dp_bits.reset(new boost::dynamic_bitset<>(size));
        // RDKit✔️✔️:   d_numOnBits = 0;
        // RDKit✔️✔️: }
        //
        // Local review: one zeroed allocation of ceil(size/64) words, exactly
        // the dynamic_bitset allocation; O(words) like the source.
        Self {
            words: vec![0; n_bits.div_ceil(64) as usize],
            n_bits,
            num_on: 0,
        }
    }

    /// Size constructor with all bits set (`bitsSet == true`).
    ///
    /// Source: `ExplicitBitVect::ExplicitBitVect(unsigned int size, bool
    /// bitsSet)` (ExplicitBitVect.cpp:23-32).
    #[must_use]
    pub fn new_filled(n_bits: u32) -> Self {
        // RDKit✔️✔️: ExplicitBitVect::ExplicitBitVect(unsigned int size, bool bitsSet) {
        // RDKit✔️✔️:   d_size = 0;
        // RDKit✔️✔️:   dp_bits = nullptr;
        // RDKit✔️✔️:   d_numOnBits = 0;
        // RDKit✔️✔️:   _initForSize(size);
        // RDKit✔️✔️:   if (bitsSet) {
        // RDKit✔️✔️:     dp_bits->set();  // set all bits to 1
        // RDKit✔️✔️:     d_numOnBits = size;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        //
        // `set()` on a dynamic_bitset only sets bits below size (tail stays
        // zero); reproduced by filling each full word and masking the tail.
        let mut fp = Self::new(n_bits);
        let full_words = (n_bits / 64) as usize;
        for word in &mut fp.words[..full_words] {
            *word = u64::MAX;
        }
        let tail = (n_bits % 64) as u32;
        if tail != 0 {
            let last = fp.words.len() - 1;
            fp.words[last] = (1u64 << tail) - 1;
        }
        fp.num_on = n_bits;
        fp
    }

    /// Constructor from explicit on-bit indices (COSMolKit extension;
    /// no DataStructs counterpart). Out-of-range
    /// indices are a structured error, not a panic.
    #[must_use]
    pub fn from_on_bits<I: IntoIterator<Item = u32>>(
        n_bits: u32,
        on_bits: I,
    ) -> Result<Self, FingerprintError> {
        let mut fp = Self::new(n_bits);
        for bit in on_bits {
            fp.set_bit(bit)?;
        }
        Ok(fp)
    }

    /// Number of bits (the length of the vector).
    ///
    /// Source: `unsigned int ExplicitBitVect::getNumBits() const { return
    /// d_size; };` (ExplicitBitVect.cpp:179).
    #[must_use]
    pub const fn n_bits(&self) -> u32 {
        // RDKit✔️✔️: unsigned int ExplicitBitVect::getNumBits() const { return d_size; };
        self.n_bits
    }

    /// Cached number of on bits.
    ///
    /// Source: `unsigned int ExplicitBitVect::getNumOnBits() const { return
    /// d_numOnBits; };` (ExplicitBitVect.cpp:180).
    #[must_use]
    pub const fn num_on_bits(&self) -> u32 {
        // RDKit✔️✔️: unsigned int ExplicitBitVect::getNumOnBits() const { return d_numOnBits; };
        self.num_on
    }

    /// Number of off bits.
    ///
    /// Source: `unsigned int ExplicitBitVect::getNumOffBits() const {
    /// return d_size - d_numOnBits; }` (ExplicitBitVect.cpp:181-183). The
    /// subtraction is unsigned `int`-width arithmetic and wraps when the
    /// (possibly stale) cache exceeds the length — reproduced with the
    /// same explicit wrapping subtraction the source performs, so the
    /// result never depends on the build profile's overflow checks.
    #[must_use]
    pub const fn num_off_bits(&self) -> u32 {
        // RDKit✔️✔️: unsigned int ExplicitBitVect::getNumOffBits() const {
        // RDKit✔️✔️:   return d_size - d_numOnBits;
        // RDKit✔️✔️: }
        off_count(self.n_bits, self.num_on)
    }

    /// Query one bit.
    ///
    /// Source: `bool ExplicitBitVect::getBit(const unsigned int which)
    /// const` (ExplicitBitVect.cpp:108-113).
    pub fn get_bit(&self, bit: u32) -> Result<bool, FingerprintError> {
        // RDKit✔️✔️: bool ExplicitBitVect::getBit(const unsigned int which) const {
        // RDKit✔️✔️:   if (which >= d_size) {
        // RDKit✔️✔️:     throw IndexErrorException(which);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return ((bool)(*dp_bits)[which]);
        // RDKit✔️✔️: }
        if bit >= self.n_bits {
            return Err(FingerprintError::SparseIndexOutOfRange {
                index: u64::from(bit),
                size: u64::from(self.n_bits),
            });
        }
        Ok(self.words[(bit / 64) as usize] & (1u64 << (bit % 64)) != 0)
    }

    /// Set one bit; returns its original state.
    ///
    /// Source: `bool ExplicitBitVect::setBit(const unsigned int which)`
    /// (ExplicitBitVect.cpp:84-95).
    pub fn set_bit(&mut self, bit: u32) -> Result<bool, FingerprintError> {
        // RDKit✔️✔️: bool ExplicitBitVect::setBit(const unsigned int which) {
        // RDKit✔️✔️:   if (which >= d_size) {
        // RDKit✔️✔️:     throw IndexErrorException(which);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if ((bool)(*dp_bits)[which]) {
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     (*dp_bits)[which] = 1;
        // RDKit✔️✔️:     ++d_numOnBits;
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        if bit >= self.n_bits {
            return Err(FingerprintError::SparseIndexOutOfRange {
                index: u64::from(bit),
                size: u64::from(self.n_bits),
            });
        }
        let mask = 1u64 << (bit % 64);
        let word = &mut self.words[(bit / 64) as usize];
        if *word & mask != 0 {
            Ok(true)
        } else {
            *word |= mask;
            self.num_on = cached_count_increment(self.num_on);
            Ok(false)
        }
    }

    /// Unset one bit; returns its original state.
    ///
    /// Source: `bool ExplicitBitVect::unsetBit(const unsigned int which)`
    /// (ExplicitBitVect.cpp:96-107).
    pub fn unset_bit(&mut self, bit: u32) -> Result<bool, FingerprintError> {
        // RDKit✔️✔️: bool ExplicitBitVect::unsetBit(const unsigned int which) {
        // RDKit✔️✔️:   if (which >= d_size) {
        // RDKit✔️✔️:     throw IndexErrorException(which);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if ((bool)(*dp_bits)[which]) {
        // RDKit✔️✔️:     (*dp_bits)[which] = 0;
        // RDKit✔️✔️:     --d_numOnBits;
        // RDKit✔️✔️:     return true;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        if bit >= self.n_bits {
            return Err(FingerprintError::SparseIndexOutOfRange {
                index: u64::from(bit),
                size: u64::from(self.n_bits),
            });
        }
        let mask = 1u64 << (bit % 64);
        let word = &mut self.words[(bit / 64) as usize];
        if *word & mask != 0 {
            *word &= !mask;
            self.num_on = cached_count_decrement(self.num_on);
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Ascending indices of all on bits.
    ///
    /// Source: `void ExplicitBitVect::getOnBits(IntVect &v) const`
    /// (ExplicitBitVect.cpp:186-197). The source clears `v`, reserves, and
    /// pushes ascending indices.
    #[must_use]
    pub fn on_bits(&self) -> Vec<u32> {
        // RDKit✔️✔️: void ExplicitBitVect::getOnBits(IntVect &v) const {
        // RDKit✔️✔️:   unsigned int nOn = getNumOnBits();
        // RDKit✔️✔️:   if (!v.empty()) {
        // RDKit✔️✔️:     IntVect().swap(v);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   v.reserve(nOn);
        // RDKit✔️✔️:   for (unsigned int i = 0; i < d_size; i++) {
        // RDKit✔️✔️:     if ((bool)(*dp_bits)[i]) {
        // RDKit✔️✔️:       v.push_back(i);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        //
        // Local review: word-wise iteration with trailing_zeros stepping is
        // asymptotically identical to the source scan (O(d_size)) with a
        // strictly smaller constant; one allocation sized to the on count,
        // like the source reserve.
        let mut out = Vec::with_capacity(self.num_on as usize);
        for (word_idx, &word) in self.words.iter().enumerate() {
            let mut remaining = word;
            while remaining != 0 {
                let offset = remaining.trailing_zeros();
                out.push((word_idx as u32) * 64 + offset);
                remaining &= remaining - 1;
            }
        }
        out
    }

    /// Clear all bits.
    ///
    /// Source: `void clearBits() override { dp_bits->reset(); }`
    /// (ExplicitBitVect.h:80). The source does **not** update `d_numOnBits`
    /// here: `getNumOnBits()` returns the stale cached value until an
    /// operation recomputes it (the binary operators and `operator~`
    /// recompute from `count()`; `setBit`/`unsetBit` increment/decrement the
    /// stale base). This upstream stale-cache semantics is reproduced
    /// exactly and pinned by regression: this is a reproduced source defect,
    /// not a divergence.
    pub fn clear_bits(&mut self) {
        // RDKit✔️✔️: void clearBits() override { dp_bits->reset(); }
        //
        // d_numOnBits is deliberately left untouched: the source assigns
        // nothing to it in clearBits, so the cached count stays stale until
        // a count-recomputing operation runs. Reproduced verbatim.
        self.words.iter_mut().for_each(|w| *w = 0);
    }

    fn same_length(&self, other: &Self) -> Result<(), FingerprintError> {
        if self.n_bits != other.n_bits {
            return Err(FingerprintError::BitLengthMismatch {
                left: u64::from(self.n_bits),
                right: u64::from(other.n_bits),
            });
        }
        Ok(())
    }

    /// Bitwise AND of two equal-length vectors.
    ///
    /// Source: `ExplicitBitVect ExplicitBitVect::operator&(const
    /// ExplicitBitVect &other) const` (ExplicitBitVect.cpp:122-127). The
    /// underlying `boost::dynamic_bitset operator&` throws
    /// `std::invalid_argument` on size mismatch; mapped to the structured
    /// `BitLengthMismatch` error for the same condition.
    pub fn and(&self, other: &Self) -> Result<Self, FingerprintError> {
        // RDKit✔️✔️: ExplicitBitVect ExplicitBitVect::operator&(const ExplicitBitVect &other) const {
        // RDKit✔️✔️:   ExplicitBitVect ans(d_size);
        // RDKit✔️✔️:   *(ans.dp_bits) = (*dp_bits) & *(other.dp_bits);
        // RDKit✔️✔️:   ans.d_numOnBits = ans.dp_bits->count();
        // RDKit✔️✔️:   return (ans);
        // RDKit✔️✔️: }
        self.same_length(other)?;
        let words = self
            .words
            .iter()
            .zip(other.words.iter())
            .map(|(a, b)| a & b);
        Ok(Self::from_words(self.n_bits, words))
    }

    /// Bitwise OR of two equal-length vectors.
    ///
    /// Source: `ExplicitBitVect ExplicitBitVect::operator|(const
    /// ExplicitBitVect &other) const` (ExplicitBitVect.cpp:129-134).
    pub fn or(&self, other: &Self) -> Result<Self, FingerprintError> {
        // RDKit✔️✔️: ExplicitBitVect ExplicitBitVect::operator|(const ExplicitBitVect &other) const {
        // RDKit✔️✔️:   ExplicitBitVect ans(d_size);
        // RDKit✔️✔️:   *(ans.dp_bits) = (*dp_bits) | *(other.dp_bits);
        // RDKit✔️✔️:   ans.d_numOnBits = ans.dp_bits->count();
        // RDKit✔️✔️:   return (ans);
        // RDKit✔️✔️: }
        self.same_length(other)?;
        let words = self
            .words
            .iter()
            .zip(other.words.iter())
            .map(|(a, b)| a | b);
        Ok(Self::from_words(self.n_bits, words))
    }

    /// Bitwise XOR of two equal-length vectors.
    ///
    /// Source: `ExplicitBitVect ExplicitBitVect::operator^(const
    /// ExplicitBitVect &other) const` (ExplicitBitVect.cpp:115-120).
    pub fn xor(&self, other: &Self) -> Result<Self, FingerprintError> {
        // RDKit✔️✔️: ExplicitBitVect ExplicitBitVect::operator^(const ExplicitBitVect &other) const {
        // RDKit✔️✔️:   ExplicitBitVect ans(d_size);
        // RDKit✔️✔️:   *(ans.dp_bits) = (*dp_bits) ^ *(other.dp_bits);
        // RDKit✔️✔️:   ans.d_numOnBits = ans.dp_bits->count();
        // RDKit✔️✔️:   return (ans);
        // RDKit✔️✔️: }
        self.same_length(other)?;
        let words = self
            .words
            .iter()
            .zip(other.words.iter())
            .map(|(a, b)| a ^ b);
        Ok(Self::from_words(self.n_bits, words))
    }

    /// Complement of the vector.
    ///
    /// Source: `ExplicitBitVect ExplicitBitVect::operator~() const`
    /// (ExplicitBitVect.cpp:154-159). `dynamic_bitset::~` flips only bits
    /// below size (tail bits stay zero).
    #[must_use]
    pub fn not(&self) -> Self {
        // RDKit✔️✔️: ExplicitBitVect ExplicitBitVect::operator~() const {
        // RDKit✔️✔️:   ExplicitBitVect ans(d_size);
        // RDKit✔️✔️:   *(ans.dp_bits) = ~(*dp_bits);
        // RDKit✔️✔️:   ans.d_numOnBits = ans.dp_bits->count();
        // RDKit✔️✔️:   return (ans);
        // RDKit✔️✔️: };
        let mut words: Vec<u64> = self.words.iter().map(|w| !w).collect();
        let tail = (self.n_bits % 64) as u32;
        if tail != 0 {
            let last = words.len() - 1;
            words[last] &= (1u64 << tail) - 1;
        }
        // Source recomputes the count from the complemented bitset
        // (`ans.d_numOnBits = ans.dp_bits->count()`), deliberately ignoring
        // a possibly stale d_numOnBits; F03 regression pins this against the
        // former `n_bits - num_on` derivation that inherited stale state.
        let num_on = words
            .iter()
            .map(|w| w.count_ones())
            .fold(0u32, |acc, c| acc + c);
        Self {
            words,
            n_bits: self.n_bits,
            num_on,
        }
    }

    /// Concatenation `self + other`.
    ///
    /// Source: `ExplicitBitVect &ExplicitBitVect::operator+=(const
    /// ExplicitBitVect &other)` (ExplicitBitVect.cpp:161-177). The size
    /// addition wraps in 32-bit `unsigned int` arithmetic; when wrapping
    /// makes a shifted target land beyond the resized length, the source's
    /// `setBit` raises `IndexErrorException`, so this method is fallible
    /// for that edge only.
    pub fn concat(&self, other: &Self) -> Result<Self, FingerprintError> {
        // RDKit✔️✔️: ExplicitBitVect &ExplicitBitVect::operator+=(const ExplicitBitVect &other) {
        // RDKit✔️✔️:   dp_bits->resize(d_size + other.d_size);
        // RDKit✔️✔️:   unsigned int original_size = d_size;
        // RDKit✔️✔️:   d_size = dp_bits->size();
        // RDKit✔️✔️:   for (unsigned i = 0; i < other.d_size; i++) {
        // RDKit✔️✔️:     if (other[i]) {
        // RDKit✔️✔️:       setBit(i + original_size);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   d_numOnBits = dp_bits->count();
        // RDKit✔️✔️:   return *this;
        // RDKit✔️✔️: }
        //
        // Local review: builds the concatenated words by shifting the first
        // partial word of `other` instead of per-bit setBit; O(words) versus
        // the source's O(other.d_size + popcount(other)) bit loop — never
        // more work per element, no extra allocations beyond the result.
        // Size arithmetic: the source adds two `unsigned int` values passed
        // to `resize` (size_t), i.e. wrapping 32-bit addition; reproduced
        // with wrapping_add (F08 correction of a previous panicking
        // checked_add that was not source behavior).
        let true_sum = u64::from(self.n_bits) + u64::from(other.n_bits);
        let new_bits = self.n_bits.wrapping_add(other.n_bits);
        let mut out = Self::new(new_bits);
        if true_sum == u64::from(new_bits) {
            // no wrap: every shifted target is provably inside new_bits, so
            // the word-shift fast path cannot hit the source's setBit bounds
            // check
            out.words[..self.words.len()].copy_from_slice(&self.words);
            let shift = (self.n_bits % 64) as u32;
            if shift == 0 {
                out.words[self.words.len()..].copy_from_slice(&other.words);
            } else {
                for (idx, &word) in other.words.iter().enumerate() {
                    let target = self.words.len() + idx;
                    out.words[target - 1] |= word << shift;
                    if target < out.words.len() {
                        out.words[target] |= word >> (64 - shift);
                    }
                }
            }
        } else {
            // 32-bit wrap: resize truncates self's bits to the wrapped
            // length and each per-bit target is computed with wrapping
            // `unsigned int` arithmetic exactly like the source loop; a
            // wrapped target beyond the new length reproduces the source's
            // setBit IndexErrorException.
            for bit in self.on_bits() {
                if bit < new_bits {
                    out.set_bit(bit)?;
                }
            }
            for bit in other.on_bits() {
                let target = self.n_bits.wrapping_add(bit);
                if target >= new_bits {
                    return Err(FingerprintError::SparseIndexOutOfRange {
                        index: u64::from(target),
                        size: u64::from(new_bits),
                    });
                }
                out.set_bit(target)?;
            }
        }
        out.num_on = out
            .words
            .iter()
            .map(|w| w.count_ones())
            .fold(0u32, |acc, c| acc + c);
        Ok(out)
    }

    fn from_words<I: Iterator<Item = u64>>(n_bits: u32, words: I) -> Self {
        let words: Vec<u64> = words.collect();
        let num_on = words
            .iter()
            .map(|w| w.count_ones())
            .fold(0u32, |a, c| a + c);
        Self {
            words,
            n_bits,
            num_on,
        }
    }
}

#[cfg(test)]
mod c2_boundary_tests {
    //! Bounded owning-module tests for the source-defined unsigned-u32
    //! arithmetic sites (C2). Boundaries are exercised directly through
    //! the production helpers — no billion-iteration sequences, no public
    //! test hooks.
    use super::{Fingerprint, cached_count_decrement, cached_count_increment, off_count};

    #[test]
    fn cached_increment_wraps_at_u32_max() {
        assert_eq!(cached_count_increment(u32::MAX), 0);
        assert_eq!(cached_count_increment(0), 1);
        assert_eq!(cached_count_increment(u32::MAX - 1), u32::MAX);
    }

    #[test]
    fn cached_decrement_wraps_at_zero() {
        assert_eq!(cached_count_decrement(0), u32::MAX);
        assert_eq!(cached_count_decrement(1), 0);
        assert_eq!(cached_count_decrement(u32::MAX), u32::MAX - 1);
    }

    #[test]
    fn off_count_wraps_when_cache_exceeds_length() {
        assert_eq!(off_count(1, 2), u32::MAX);
        assert_eq!(off_count(0, 1), u32::MAX);
        assert_eq!(off_count(70, 6), 64);
        assert_eq!(off_count(u32::MAX, 0), u32::MAX);
    }

    #[test]
    fn public_stale_sequence_matches_supervisor_probe() {
        // bounded end-to-end shape: new(1), set, clear, set -> cache 2,
        // off u32::MAX, equality with fresh one-bit vector restored by C1
        let mut one = Fingerprint::new(1);
        one.set_bit(0).unwrap();
        one.clear_bits();
        one.set_bit(0).unwrap();
        assert_eq!(one.num_on_bits(), 2);
        assert_eq!(one.num_off_bits(), u32::MAX);
        assert_eq!(one, Fingerprint::new(1).with_bit_zero());
    }

    impl Fingerprint {
        /// test-only tiny constructor: fresh vector with bit 0 set and a
        /// truthful cache
        fn with_bit_zero(mut self) -> Self {
            self.set_bit(0).unwrap();
            self
        }
    }
}
