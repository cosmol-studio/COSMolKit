//! Scalar hash primitives used by fingerprint values: source-backed port
//! of the pinned RDKit `Code/RDGeneral/hash` integral closure
//! (`hash_fwd.hpp`, `hash.hpp`).
//!
//! `std::hash_result_t` is `std::uint32_t` (hash_fwd.hpp:18-20): every
//! seed, intermediate and result below is 32-bit with wrapping arithmetic,
//! exactly as declared by the source. Under the pinned ABI (gcc x86-64),
//! `long`/`unsigned long` are 64-bit and use the **plain cast** overloads
//! (low-32 truncation), distinct from the 64-bit mixing helpers used by
//! `long long`/`unsigned long long` (F22/F23). No `DefaultHasher` or
//! replacement hash algorithm appears anywhere in this module.

/// `hash_result_t` (hash_fwd.hpp:18-20).
pub type HashResult = u32;

/// `hash_value(bool)` (hash.hpp:138-140).
#[must_use]
pub const fn hash_value_bool(v: bool) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(bool v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v as HashResult
}

/// `hash_value(char/signed char/short/int)` signed small-type overload:
/// value-converted to `int` then cast; negative values occupy the full
/// two's-complement 32-bit pattern (hash.hpp:142-164).
#[must_use]
pub const fn hash_value_i8(v: i8) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(signed char v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v as HashResult
}

/// `hash_value(unsigned char)` (hash.hpp:146-148).
#[must_use]
pub const fn hash_value_u8(v: u8) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(unsigned char v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v as HashResult
}

/// `hash_value(short)` (hash.hpp:154-156).
#[must_use]
pub const fn hash_value_i16(v: i16) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(short v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v as HashResult
}

/// `hash_value(unsigned short)` (hash.hpp:158-160).
#[must_use]
pub const fn hash_value_u16(v: u16) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(unsigned short v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v as HashResult
}

/// `hash_value(int)` (hash.hpp:162-164): two's-complement cast.
#[must_use]
pub const fn hash_value_i32(v: i32) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(int v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v as HashResult
}

/// `hash_value(unsigned int)` (hash.hpp:166-168): identity.
#[must_use]
pub const fn hash_value_u32(v: u32) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(unsigned int v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v
}

/// `hash_detail::hash_value_unsigned<T>` instantiated for
/// `boost::ulong_long_type` (64-bit), the overload reached by
/// `hash_value(unsigned long long)` (hash.hpp:118-135, 189-191).
///
/// For 64-bit `T`, `length = (64-1)/32 = 1`, so the loop folds the high
/// 32 bits once and then the low 32 bits, with all arithmetic wrapping in
/// the declared unsigned `hash_result_t`. Distinct from the plain-cast
/// `unsigned long` overload (F21): `0x1_0000_0000u64` mixes to `65` here
/// but truncates to `0` there.
#[must_use]
pub const fn hash_value_u64(val: u64) -> HashResult {
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: inline std::hash_result_t hash_value_unsigned(T val) {
    // RDKit✔️✔️:   const int hash_result_t_bits =
    // RDKit✔️✔️:       std::numeric_limits<std::hash_result_t>::digits;
    // RDKit✔️✔️:   // ceiling(std::numeric_limits<T>::digits / hash_result_t_bits) - 1
    // RDKit✔️✔️:   const int length = (std::numeric_limits<T>::digits - 1) / hash_result_t_bits;
    // RDKit✔️✔️:   std::hash_result_t seed = 0;
    // RDKit✔️✔️:   // Hopefully, this loop can be unrolled.
    // RDKit✔️✔️:   for (unsigned int i = length * hash_result_t_bits; i > 0;
    // RDKit✔️✔️:        i -= hash_result_t_bits) {
    // RDKit✔️✔️:     seed ^= (std::hash_result_t)(val >> i) + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   seed ^= (std::hash_result_t)val + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️:   return seed;
    // RDKit✔️✔️: }
    let mut seed: HashResult = 0;
    let mut i = 32u32;
    while i > 0 {
        seed ^= ((val >> i) as HashResult)
            .wrapping_add(seed << 6)
            .wrapping_add(seed >> 2);
        i -= 32;
    }
    seed ^= (val as HashResult)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2);
    seed
}

/// `hash_value(unsigned long)` (hash.hpp:174-176). Under the pinned ABI
/// `unsigned long` is 64-bit, and this overload is the **plain cast**
/// (defined 64→32 truncation), unlike `unsigned long long` which uses the
/// mixing helper (F22). Pinned distinction: `0x1_0000_0001usize` hashes to
/// `1` here while the u64 overload mixes both halves.
#[must_use]
pub const fn hash_value_usize(v: usize) -> HashResult {
    // RDKit✔️✔️: inline std::hash_result_t hash_value(unsigned long v) {
    // RDKit✔️✔️:   return static_cast<std::hash_result_t>(v);
    // RDKit✔️✔️: }
    v as HashResult
}

/// `hash_detail::hash_value_signed<T>` instantiated for
/// `boost::long_long_type` (64-bit), the overload reached by
/// `hash_value(long long)` (hash.hpp:98-116, 185-187).
///
/// For 64-bit `T`, `length = (63-1)/32 = 1`: the loop folds
/// `positive >> 32` once (`positive = val < 0 ? -1 - val : val`, which is
/// always representable), then folds `(hash_result_t)val` — the
/// two's-complement low word of the original signed value. All arithmetic
/// wraps in the declared unsigned `hash_result_t`.
#[must_use]
pub const fn hash_value_i64(val: i64) -> HashResult {
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: inline std::hash_result_t hash_value_signed(T val) {
    // RDKit✔️✔️:   const int hash_result_t_bits =
    // RDKit✔️✔️:       std::numeric_limits<std::hash_result_t>::digits;
    // RDKit✔️✔️:   // ceiling(std::numeric_limits<T>::digits / hash_result_t_bits) - 1
    // RDKit✔️✔️:   const int length = (std::numeric_limits<T>::digits - 1) / hash_result_t_bits;
    // RDKit✔️✔️:   std::hash_result_t seed = 0;
    // RDKit✔️✔️:   T positive = val < 0 ? -1 - val : val;
    // RDKit✔️✔️:   // Hopefully, this loop can be unrolled.
    // RDKit✔️✔️:   for (unsigned int i = length * hash_result_t_bits; i > 0;
    // RDKit✔️✔️:        i -= hash_result_t_bits) {
    // RDKit✔️✔️:     seed ^= (std::hash_result_t)(positive >> i) + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   seed ^= (std::hash_result_t)val + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️:   return seed;
    // RDKit✔️✔️: }
    let positive: u64 = if val < 0 {
        // `-1 - val` is representable for every negative i64 (max is i64::MAX)
        (-1 - val) as u64
    } else {
        val as u64
    };
    let mut seed: HashResult = 0;
    let mut i = 32u32;
    while i > 0 {
        seed ^= ((positive >> i) as HashResult)
            .wrapping_add(seed << 6)
            .wrapping_add(seed >> 2);
        i -= 32;
    }
    seed ^= (val as HashResult)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2);
    seed
}

/// `gboost::hash_combine` for the fingerprint-value call-site overload
/// (`T = unsigned int`, whose `hash_value` is the identity), from
/// hash.hpp:209-219. The seed and result are `std::hash_result_t`
/// (`uint32_t`, hash_fwd.hpp:18-20); all arithmetic wraps in that
/// declared unsigned width.
pub fn hash_combine(seed: &mut HashResult, value: u32) {
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: inline void hash_combine(std::hash_result_t& seed, T const& v)
    // RDKit✔️✔️: {
    // RDKit✔️✔️:   gboost::hash<T> hasher;
    // RDKit✔️✔️:   seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    // RDKit✔️✔️: }
    *seed ^= value
        .wrapping_add(0x9e3779b9)
        .wrapping_add(seed.wrapping_shl(6))
        .wrapping_add(seed.wrapping_shr(2));
}

/// `hash_range(first, last)` for the unsigned-int value call sites:
/// ordered accumulation from an initial zero seed (hash.hpp:221-230).
#[must_use]
pub fn hash_range(values: &[u32]) -> HashResult {
    // RDKit✔️✔️: template <typename It>
    // RDKit✔️✔️: inline std::hash_result_t hash_range(It first, It last) {
    // RDKit✔️✔️:   std::hash_result_t seed = 0;
    // RDKit✔️✔️:   for (; first != last; ++first) {
    // RDKit✔️✔️:     hash_combine(seed, *first);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return seed;
    // RDKit✔️✔️: }
    let mut seed = 0;
    hash_range_seeded(&mut seed, values);
    seed
}

/// `hash_range(seed&, first, last)`: ordered accumulation into a
/// caller-supplied seed (hash.hpp:232-237).
pub fn hash_range_seeded(seed: &mut HashResult, values: &[u32]) {
    // RDKit✔️✔️: template <typename It>
    // RDKit✔️✔️: inline void hash_range(std::hash_result_t& seed, It first, It last) {
    // RDKit✔️✔️:   for (; first != last; ++first) {
    // RDKit✔️✔️:     hash_combine(seed, *first);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for &value in values {
        hash_combine(seed, value);
    }
}
