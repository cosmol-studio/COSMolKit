//! Sparse count vectors: source-backed port of RDKit
//! `Code/DataStructs/SparseIntVect.h` in both index specializations used
//! by the fingerprint value boundary — `SparseIntVect<std::uint64_t>`
//! (`SparseCountFingerprint`) and `SparseIntVect<std::uint32_t>`
//! (`SparseCountFingerprint32`), each with signed `int` (i32) values
//! preserved exactly.
//!
//! Signed arithmetic is checked only at the supervisor-approved executed
//! UB sites (overflow, `abs(MIN)`, executed division by zero,
//! `MIN / -1`, right-only negation of MIN), reported as structured
//! `UndefinedArithmetic`; every defined source result is exact.

use std::collections::BTreeMap;

use crate::FingerprintError;

macro_rules! sparse_count_storage {
    ($ty:ident, $width:ty) => {
        impl $ty {
            /// Length constructor.
            ///
            /// Source: `SparseIntVect(IndexType length) : d_length(length)
            /// {}` (SparseIntVect.h:34-35).
            #[must_use]
            pub const fn new(length: $width) -> Self {
                // RDKit✔️✔️: SparseIntVect(IndexType length) : d_length(length) {}
                Self {
                    length,
                    data: BTreeMap::new(),
                }
            }

            /// Vector length.
            ///
            /// Source: `IndexType getLength() const { return d_length; }`
            /// (SparseIntVect.h:110).
            #[must_use]
            pub const fn length(&self) -> $width {
                // RDKit✔️✔️: IndexType getLength() const { return d_length; }
                self.length
            }

            /// Source index check (per-width `numeric_limits<IndexType>`).
            ///
            /// Source: `bool checkIndex(IndexType idx) const { if (idx < 0
            /// || idx > d_length || (idx == d_length && d_length <
            /// std::numeric_limits<IndexType>::max())) { return false; }
            /// return true; }` (SparseIntVect.h:413-419). For the unsigned
            /// instantiations `idx < 0` is never true.
            const fn check_index(&self, idx: $width) -> bool {
                // RDKit✔️✔️: bool checkIndex(IndexType idx) const {
                // RDKit✔️✔️:   if (idx < 0 || idx > d_length ||
                // RDKit✔️✔️:       (idx == d_length && d_length < std::numeric_limits<IndexType>::max())) {
                // RDKit✔️✔️:     return false;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return true;
                // RDKit✔️✔️: }
                !(idx > self.length || (idx == self.length && self.length < <$width>::MAX))
            }

            fn checked(&self, idx: $width) -> Result<(), FingerprintError> {
                if !self.check_index(idx) {
                    return Err(FingerprintError::SparseIndexOutOfRange {
                        index: u64::from(idx),
                        size: u64::from(self.length),
                    });
                }
                Ok(())
            }

            /// Value at an index (0 when absent).
            ///
            /// Source: `int getVal(IndexType idx) const`
            /// (SparseIntVect.h:77-87).
            pub fn value(&self, idx: $width) -> Result<i32, FingerprintError> {
                // RDKit✔️✔️: int getVal(IndexType idx) const {
                // RDKit✔️✔️:   if (!checkIndex(idx)) {
                // RDKit✔️✔️:     throw IndexErrorException(static_cast<int>(idx));
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   int res = 0;
                // RDKit✔️✔️:   typename StorageType::const_iterator iter = d_data.find(idx);
                // RDKit✔️✔️:   if (iter != d_data.end()) {
                // RDKit✔️✔️:     res = iter->second;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return res;
                // RDKit✔️✔️: }
                self.checked(idx)?;
                Ok(self.data.get(&idx).copied().unwrap_or(0))
            }

            /// Set the value at an index; zero erases the entry.
            ///
            /// Source: `void setVal(IndexType idx, int val)`
            /// (SparseIntVect.h:90-99).
            pub fn set_value(&mut self, idx: $width, val: i32) -> Result<(), FingerprintError> {
                // RDKit✔️✔️: void setVal(IndexType idx, int val) {
                // RDKit✔️✔️:   if (!checkIndex(idx)) {
                // RDKit✔️✔️:     throw IndexErrorException(static_cast<int>(idx));
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   if (val != 0) {
                // RDKit✔️✔️:     d_data[idx] = val;
                // RDKit✔️✔️:   } else {
                // RDKit✔️✔️:     d_data.erase(idx);
                // RDKit✔️✔️:   }
                // RDKit✔️✔️: }
                self.checked(idx)?;
                if val != 0 {
                    self.data.insert(idx, val);
                } else {
                    self.data.remove(&idx);
                }
                Ok(())
            }

            /// Nonzero entries in ascending index order.
            ///
            /// Source: `const StorageType &getNonzeroElements() const {
            /// return d_data; }` (SparseIntVect.h:130). Entries whose
            /// stored value is zero (producible only through the scalar
            /// operators, which do not erase) remain visible, exactly like
            /// the source map.
            #[must_use]
            pub const fn nonzero_elements(&self) -> &BTreeMap<$width, i32> {
                // RDKit✔️✔️: const StorageType &getNonzeroElements() const { return d_data; }
                &self.data
            }

            fn same_length(&self, other: &Self) -> Result<(), FingerprintError> {
                if other.length != self.length {
                    return Err(FingerprintError::BitLengthMismatch {
                        left: u64::from(self.length),
                        right: u64::from(other.length),
                    });
                }
                Ok(())
            }

            /// Fuzzy intersection: per-index minimum.
            ///
            /// Source: `SparseIntVect &operator&=(const SparseIntVect
            /// &other)` (SparseIntVect.h:135-164). Own entries absent in
            /// `other` are dropped; common keys keep the minimum.
            ///
            /// One result-map clone followed by ordered, cursor-based retention.
            /// Aggregate O(n + m), including tree repair and internal-key
            /// deletion; see performance_alignment.md section 2 for the
            /// height-from-leaf event bound. Result storage is O(n), without
            /// a second input-sized buffer. Heap allocation counts are not
            /// identical to std::map: their node capacities differ.
            pub fn fuzzy_and(&self, other: &Self) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator&=(const SparseIntVect<IndexType> &other) {
                // RDKit✔️✔️:   if (other.d_length != d_length) {
                // RDKit✔️✔️:     throw ValueErrorException("SparseIntVect size mismatch");
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   typename StorageType::const_iterator oIter = other.d_data.begin();
                // RDKit✔️✔️:   while (iter != d_data.end()) {
                // RDKit✔️✔️:     // we're relying on the fact that the maps are sorted:
                // RDKit✔️✔️:     while (oIter != other.d_data.end() && oIter->first < iter->first) {
                // RDKit✔️✔️:       ++oIter;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     if (oIter != other.d_data.end() && oIter->first == iter->first) {
                // RDKit✔️✔️:       // found it:
                // RDKit✔️✔️:       if (oIter->second < iter->second) {
                // RDKit✔️✔️:         iter->second = oIter->second;
                // RDKit✔️✔️:       }
                // RDKit✔️✔️:       ++oIter;
                // RDKit✔️✔️:       ++iter;
                // RDKit✔️✔️:     } else {
                // RDKit✔️✔️:       // not there; our value is zero, which means
                // RDKit✔️✔️:       // we should remove this value:
                // RDKit✔️✔️:       typename StorageType::iterator tmpIter = iter;
                // RDKit✔️✔️:       ++tmpIter;
                // RDKit✔️✔️:       d_data.erase(iter);
                // RDKit✔️✔️:       iter = tmpIter;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                // RDKit✔️✔️: const SparseIntVect<IndexType> operator&(
                // RDKit✔️✔️:     const SparseIntVect<IndexType> &other) const {
                // RDKit✔️✔️:   SparseIntVect<IndexType> res(*this);
                // RDKit✔️✔️:   return res &= other;
                // RDKit✔️✔️: }
                //
                // Clone once, then retain shared keys with the signed minimum.
                // The right iterator advances monotonically. The installed
                // retain implementation allocates no new nodes; child merges
                // consume initial nodes. At height j above leaves, separator
                // arrivals are bounded by initial slots and repairs at j/j-1.
                // The weighted geometric sum of navigation costs is O(n),
                // not O(n log n), even with predecessor replacement: the
                // replacement predecessor was already visited; traversal
                // resumes after it rather than counting it as a new arrival.
                // No claim of identical allocator counts or wall-clock times.
                self.same_length(other)?;
                let mut data = self.data.clone();
                let mut other_iter = other.data.iter();
                let mut other_cur = other_iter.next();
                data.retain(|key, own_val| {
                    while let Some((&other_key, _)) = other_cur {
                        if other_key < *key {
                            other_cur = other_iter.next();
                        } else {
                            break;
                        }
                    }
                    if let Some((&other_key, &other_val)) = other_cur {
                        if other_key == *key {
                            if other_val < *own_val {
                                *own_val = other_val;
                            }
                            other_cur = other_iter.next();
                            return true;
                        }
                    }
                    false
                });
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Per-entry scalar subtract. Zero results persist (no erase in
            /// source).
            ///
            /// Source: `SparseIntVect &operator-=(int v)` (SparseIntVect.h:312-319).
            ///
            /// Boundary (Supervisor resolution 2026-09-24): signed subtraction
            /// overflow is C++-undefined when it executes; the Rust port performs
            /// `checked_sub` on existing entries only and returns
            /// `UndefinedArithmetic` exactly then. A COSMolKit safety boundary,
            /// not an RDKit exception or all-input parity claim. Zero results are
            /// retained, absent indices stay absent, and the input is unchanged on
            /// success and error.
            pub fn with_subtracted_scalar(&self, v: i32) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator-=(int v) {
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   while (iter != d_data.end()) {
                // RDKit✔️✔️:     iter->second -= v;
                // RDKit✔️✔️:     ++iter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                //
                // Local review: clones the map once, mutates stored values in place.
                let mut data = self.data.clone();
                for value in data.values_mut() {
                    *value = value
                        .checked_sub(v)
                        .ok_or(FingerprintError::UndefinedArithmetic {
                            site: "SparseIntVect::operator-=(int)",
                        })?;
                }
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Per-entry scalar add. Zero results persist (no erase in source).
            ///
            /// Source: `SparseIntVect &operator+=(int v)` (SparseIntVect.h:300-307).
            ///
            /// Boundary (Supervisor resolution 2026-09-24): signed add overflow is
            /// C++-undefined when it executes; the Rust port performs `checked_add`
            /// on existing entries only and returns `UndefinedArithmetic` exactly
            /// then. A COSMolKit safety boundary, not an RDKit exception or
            /// all-input parity claim. Zero results are retained, absent indices
            /// stay absent (no densification), and the input is unchanged on
            /// success and error.
            pub fn with_added_scalar(&self, v: i32) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator+=(int v) {
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   while (iter != d_data.end()) {
                // RDKit✔️✔️:     iter->second += v;
                // RDKit✔️✔️:     ++iter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                //
                // Local review: clones the map once, mutates stored values in place.
                let mut data = self.data.clone();
                for value in data.values_mut() {
                    *value = value
                        .checked_add(v)
                        .ok_or(FingerprintError::UndefinedArithmetic {
                            site: "SparseIntVect::operator+=(int)",
                        })?;
                }
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Per-entry scalar divide (truncation toward zero). The source does
            /// not erase zero results.
            ///
            /// Source: `SparseIntVect &operator/=(int v)` (SparseIntVect.h:288-295).
            ///
            /// Boundary (Supervisor resolution 2026-09-24): division by zero and
            /// `i32::MIN / -1` are C++-undefined **when a division actually
            /// executes**; the Rust port performs `checked_div` and returns
            /// `UndefinedArithmetic` exactly for those executed operations. Empty
            /// storage divided by zero therefore succeeds (no division executes),
            /// while a stored entry (zero or nonzero) divided by zero errors. This
            /// is a COSMolKit safety boundary, not an RDKit exception and not a
            /// claim of all-input parity. Defined quotients truncate toward zero
            /// exactly like the source; entry set, ordering, length and the input
            /// value are unchanged on success and error.
            pub fn with_divided_scalar(&self, v: i32) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator/=(int v) {
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   while (iter != d_data.end()) {
                // RDKit✔️✔️:     iter->second /= v;
                // RDKit✔️✔️:     ++iter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                //
                // Local review: clones the map once and mutates stored values in
                // place, mirroring the source's in-place entry update.
                let mut data = self.data.clone();
                for value in data.values_mut() {
                    *value = value
                        .checked_div(v)
                        .ok_or(FingerprintError::UndefinedArithmetic {
                            site: "SparseIntVect::operator/=(int)",
                        })?;
                }
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Per-entry scalar multiply.
            ///
            /// Source: `SparseIntVect &operator*=(int v)` (SparseIntVect.h:276-283).
            /// The source does not erase entries that become zero, so multiplied-
            /// through-zero entries persist.
            ///
            /// Boundary (Supervisor resolution 2026-09-24): signed multiply overflow
            /// (including `i32::MIN * -1`) is C++-undefined; the Rust port performs
            /// `checked_mul` and returns `UndefinedArithmetic` exactly when such an
            /// operation executes. This is a COSMolKit safety boundary, not an
            /// RDKit exception and not a claim of all-input parity. Every
            /// defined-domain product is exact; entry set, ordering, length and the
            /// input value are unchanged on success and error.
            pub fn with_multiplied_scalar(&self, v: i32) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator*=(int v) {
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   while (iter != d_data.end()) {
                // RDKit✔️✔️:     iter->second *= v;
                // RDKit✔️✔️:     ++iter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                //
                // Local review (acceptance correction): the installed Rust
                // alloc BTreeMap::clone recursively clones each
                // integer-key/value node exactly once
                // (alloc::collections::btree::map::Clone::clone_subtree),
                // then values_mut walks those n nodes once — O(n) time and
                // O(n) result storage for the immutable value-style result.
                // The C++ source mutates *this in place (no allocation);
                // the Rust form is the copy-then-scalar shape required by
                // the value operator returning a new object, not the
                // source's in-place update.
                let mut data = self.data.clone();
                for value in data.values_mut() {
                    *value = value
                        .checked_mul(v)
                        .ok_or(FingerprintError::UndefinedArithmetic {
                            site: "SparseIntVect::operator*=(int)",
                        })?;
                }
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Per-index difference; entries that cancel to zero are erased and
            /// absent own entries become the negated other value.
            ///
            /// Source: `SparseIntVect &operator-=(const SparseIntVect &other)`
            /// (SparseIntVect.h:243-270).
            pub fn with_subtracted(&self, other: &Self) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator-=(const SparseIntVect<IndexType> &other) {
                // RDKit✔️✔️:   if (other.d_length != d_length) {
                // RDKit✔️✔️:     throw ValueErrorException("SparseIntVect size mismatch");
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   typename StorageType::const_iterator oIter = other.d_data.begin();
                // RDKit✔️✔️:   while (oIter != other.d_data.end()) {
                // RDKit✔️✔️:     while (iter != d_data.end() && iter->first < oIter->first) {
                // RDKit✔️✔️:       ++iter;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     if (iter != d_data.end() && oIter->first == iter->first) {
                // RDKit✔️✔️:       // found it:
                // RDKit✔️✔️:       iter->second -= oIter->second;
                // RDKit✔️✔️:       if (!iter->second) {
                // RDKit✔️✔️:         typename StorageType::iterator tIter = iter;
                // RDKit✔️✔️:         ++tIter;
                // RDKit✔️✔️:         d_data.erase(iter);
                // RDKit✔️✔️:         iter = tIter;
                // RDKit✔️✔️:       } else {
                // RDKit✔️✔️:         ++iter;
                // RDKit✔️✔️:       }
                // RDKit✔️✔️:     } else {
                // RDKit✔️✔️:       d_data[oIter->first] = -oIter->second;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     ++oIter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                // Signed `-=` overflow and `-other` negation of i32::MIN are C++
                // UB when executed: reported via `UndefinedArithmetic` per the
                // approved safety boundary (2026-09-24);
                // never silently wrapped.
                self.same_length(other)?;
                let mut data = self.data.clone();
                for (idx, &other_val) in &other.data {
                    match data.get_mut(idx) {
                        Some(own_val) => {
                            *own_val = own_val.checked_sub(other_val).ok_or(
                                FingerprintError::UndefinedArithmetic {
                                    site: "SparseIntVect::operator-= accumulation",
                                },
                            )?;
                            if *own_val == 0 {
                                data.remove(idx);
                            }
                        }
                        None => {
                            let negated = other_val.checked_neg().ok_or(
                                FingerprintError::UndefinedArithmetic {
                                    site: "SparseIntVect::operator-= negation",
                                },
                            )?;
                            data.insert(*idx, negated);
                        }
                    }
                }
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Per-index sum; entries that cancel to zero are erased.
            ///
            /// Source: `SparseIntVect &operator+=(const SparseIntVect &other)`
            /// (SparseIntVect.h:209-236).
            pub fn with_added(&self, other: &Self) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator+=(const SparseIntVect<IndexType> &other) {
                // RDKit✔️✔️:   if (other.d_length != d_length) {
                // RDKit✔️✔️:     throw ValueErrorException("SparseIntVect size mismatch");
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   typename StorageType::const_iterator oIter = other.d_data.begin();
                // RDKit✔️✔️:   while (oIter != other.d_data.end()) {
                // RDKit✔️✔️:     while (iter != d_data.end() && iter->first < oIter->first) {
                // RDKit✔️✔️:       ++iter;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     if (iter != d_data.end() && oIter->first == iter->first) {
                // RDKit✔️✔️:       // found it:
                // RDKit✔️✔️:       iter->second += oIter->second;
                // RDKit✔️✔️:       if (!iter->second) {
                // RDKit✔️✔️:         typename StorageType::iterator tIter = iter;
                // RDKit✔️✔️:         ++tIter;
                // RDKit✔️✔️:         d_data.erase(iter);
                // RDKit✔️✔️:         iter = tIter;
                // RDKit✔️✔️:       } else {
                // RDKit✔️✔️:         ++iter;
                // RDKit✔️✔️:       }
                // RDKit✔️✔️:     } else {
                // RDKit✔️✔️:       d_data[oIter->first] = oIter->second;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     ++oIter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                // Signed `+=` overflow is C++ UB when executed: reported via
                // `UndefinedArithmetic` per the resolved supervisor decision
                // (2026-09-24); never silently wrapped.
                self.same_length(other)?;
                let mut data = self.data.clone();
                for (idx, &other_val) in &other.data {
                    match data.get_mut(idx) {
                        Some(own_val) => {
                            *own_val = own_val.checked_add(other_val).ok_or(
                                FingerprintError::UndefinedArithmetic {
                                    site: "SparseIntVect::operator+= accumulation",
                                },
                            )?;
                            if *own_val == 0 {
                                data.remove(idx);
                            }
                        }
                        None => {
                            data.insert(*idx, other_val);
                        }
                    }
                }
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Fuzzy union: per-index maximum.
            ///
            /// Source: `SparseIntVect &operator|=(const SparseIntVect &other)`
            /// (SparseIntVect.h:174-202).
            ///
            /// Performance Step 10/20: one O(n) clone of the left map (the
            /// copy the pinned source's non-mutating value operator pays),
            /// then two monotonic scans — shared keys updated to the
            /// source maximum in place (zero lookups), and right-only rows
            /// (r) inserted via ordered insert. Total O(n + m + r log(n+r)),
            /// O(n+m) for fully overlapping supports including the clone;
            /// the same primitives and class as the source value operator
            /// (proof: performance_alignment.md §3; comparable
            /// order/peak/result-storage, not identical allocation
            /// counts — node capacities differ: up to 11 KVs per
            /// BTreeMap node vs 1 per std::map node).
            pub fn fuzzy_or(&self, other: &Self) -> Result<Self, FingerprintError> {
                // RDKit✔️✔️: SparseIntVect<IndexType> &operator|=(const SparseIntVect<IndexType> &other) {
                // RDKit✔️✔️:   if (other.d_length != d_length) {
                // RDKit✔️✔️:     throw ValueErrorException("SparseIntVect size mismatch");
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   typename StorageType::iterator iter = d_data.begin();
                // RDKit✔️✔️:   typename StorageType::const_iterator oIter = other.d_data.begin();
                // RDKit✔️✔️:   while (iter != d_data.end()) {
                // RDKit✔️✔️:     // we're relying on the fact that the maps are sorted:
                // RDKit✔️✔️:     while (oIter != other.d_data.end() && oIter->first < iter->first) {
                // RDKit✔️✔️:       d_data[oIter->first] = oIter->second;
                // RDKit✔️✔️:       ++oIter;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     if (oIter != other.d_data.end() && oIter->first == iter->first) {
                // RDKit✔️✔️:       // found it:
                // RDKit✔️✔️:       if (oIter->second > iter->second) {
                // RDKit✔️✔️:         iter->second = oIter->second;
                // RDKit✔️✔️:       }
                // RDKit✔️✔️:       ++oIter;
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:     ++iter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   // finish up the other vect:
                // RDKit✔️✔️:   while (oIter != other.d_data.end()) {
                // RDKit✔️✔️:     d_data[oIter->first] = oIter->second;
                // RDKit✔️✔️:     ++oIter;
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return *this;
                // RDKit✔️✔️: }
                //
                // RDKit✔️✔️: const SparseIntVect<IndexType> operator|(
                // RDKit✔️✔️:     const SparseIntVect<IndexType> &other) const {
                // RDKit✔️✔️:   SparseIntVect<IndexType> res(*this);
                // RDKit✔️✔️:   return res |= other;
                // RDKit✔️✔️: }
                //
                // Performance Step 10 (2026-09-24): the value-style result mirrors
                // the source's `const operator|` shape — one O(n) clone of the
                // left map, then two monotonic scans. The first pairs mutable
                // clone entries with right entries and updates shared keys to
                // the source maximum in place. The second walks the ORIGINAL
                // immutable left keys and the right entries side by side and
                // inserts ONLY right-exclusive rows into the result (r ordered
                // insertions, O(r log(n+r))). Iterator state is O(1); there is
                // no get/get_mut/contains_key per shared key, no per-key
                // search of the original left, no full result rebuild, and no
                // input-sized buffer beyond the result. Aggregate:
                // O(n + m + r log(n+r)); for fully overlapping n-key supports
                // (r=0) this is O(n + m) including the clone.
                // This is source-based complexity review, not timing evidence.
                self.same_length(other)?;
                let mut data = self.data.clone();
                // Pass 1: monotonic shared-key maximum updates in place.
                let mut other_iter = other.data.iter();
                let mut other_cur = other_iter.next();
                for (_key, own_val) in data.iter_mut() {
                    while let Some((&other_key, _)) = other_cur {
                        if other_key < *_key {
                            other_cur = other_iter.next();
                        } else {
                            break;
                        }
                    }
                    if let Some((&other_key, &other_val)) = other_cur {
                        if other_key == *_key {
                            if other_val > *own_val {
                                *own_val = other_val;
                            }
                            other_cur = other_iter.next();
                        }
                    }
                }
                // Pass 2: monotonic scan over the ORIGINAL left keys and the
                // right entries, inserting only right-exclusive rows. Because
                // both sequences are ordered and distinct-keyed, a right key
                // equals a left key iff it is the one the paired walk is
                // currently facing; otherwise the ordering decides skip vs
                // insert.
                let mut left_iter = self.data.iter();
                let mut left_cur = left_iter.next();
                for (&other_key, &other_val) in &other.data {
                    while let Some((&left_key, _)) = left_cur {
                        if left_key < other_key {
                            left_cur = left_iter.next();
                        } else {
                            break;
                        }
                    }
                    match left_cur {
                        Some((&left_key, _)) if left_key == other_key => {
                            // shared key already updated in pass 1
                            left_cur = left_iter.next();
                        }
                        _ => {
                            data.insert(other_key, other_val);
                        }
                    }
                }
                Ok(Self {
                    length: self.length,
                    data,
                })
            }

            /// Sum of all elements.
            ///
            /// Source: `int getTotalVal(bool doAbs = false) const`
            /// (SparseIntVect.h:114-125). The accumulator is the source's `int`
            /// width. Signed accumulation overflow and `abs(i32::MIN)` are C++ UB
            /// when executed; under the approved safety boundary (2026-09-24)
            /// those executed edges are **reported**
            /// (`UndefinedArithmetic`) — a CK safety boundary, not an RDKit
            /// exception — and every defined accumulation is exact.
            pub fn total_value(&self, do_abs: bool) -> Result<i32, FingerprintError> {
                // RDKit✔️✔️: int getTotalVal(bool doAbs = false) const {
                // RDKit✔️✔️:   int res = 0;
                // RDKit✔️✔️:   typename StorageType::const_iterator iter;
                // RDKit✔️✔️:   for (iter = d_data.begin(); iter != d_data.end(); ++iter) {
                // RDKit✔️✔️:     if (!doAbs) {
                // RDKit✔️✔️:       res += iter->second;
                // RDKit✔️✔️:     } else {
                // RDKit✔️✔️:       res += abs(iter->second);
                // RDKit✔️✔️:     }
                // RDKit✔️✔️:   }
                // RDKit✔️✔️:   return res;
                // RDKit✔️✔️: }
                let mut res: i32 = 0;
                for &value in self.data.values() {
                    let addend = if do_abs {
                        value
                            .checked_abs()
                            .ok_or(FingerprintError::UndefinedArithmetic {
                                site: "SparseIntVect::getTotalVal abs",
                            })?
                    } else {
                        value
                    };
                    res = res
                        .checked_add(addend)
                        .ok_or(FingerprintError::UndefinedArithmetic {
                            site: "SparseIntVect::getTotalVal accumulation",
                        })?;
                }
                Ok(res)
            }
        }
    };
}

/// Sparse integer-count vector (`SparseIntVect<std::uint64_t>`
/// equivalent).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseCountFingerprint {
    length: u64,
    data: BTreeMap<u64, i32>,
}

/// Sparse integer-count vector (`SparseIntVect<std::uint32_t>`
/// equivalent). Same source template, 32-bit index width: the maximum
/// valid index is `u32::MAX`, and `checkIndex` permits `idx == length`
/// only when `length == u32::MAX`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SparseCountFingerprint32 {
    length: u32,
    data: BTreeMap<u32, i32>,
}

sparse_count_storage!(SparseCountFingerprint, u64);
sparse_count_storage!(SparseCountFingerprint32, u32);
