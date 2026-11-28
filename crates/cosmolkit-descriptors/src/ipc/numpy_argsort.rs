//! NumPy2.5.3 float64 ArgQSort, fixed X86_V4/AVX512 eight-lane order.
//! Intel x86-simd-sort fa944efbea1a33426b3f8a9a21e23794ab8aa300 (BSD-3-Clause).
//! Copyright (C) 2022 Intel Corporation; source notice retained in LICENSE.
//! Scalar lane storage preserves the named reference's equivalent-key ordering
//! on native and WASM; this is not a Rust stable/unstable-sort substitution.
//! Scalar register copies/comparisons preserve the SIMD dependency graph but
//! have known extra constant cost relative to the reference AVX512 kernels.

#[derive(Clone, Copy)]
struct Lane {
    key: f64,
    index: usize,
}
const EMPTY: Lane = Lane {
    key: f64::INFINITY,
    index: usize::MAX,
};
type Register = [Lane; 8];

fn compare_merge(values: &mut Register, other: Register, mask: u8) {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename vtype1,
    // NumPy✔️❌:           typename vtype2,
    // NumPy✔️❌:           typename reg_t1 = typename vtype1::reg_t,
    // NumPy✔️❌:           typename reg_t2 = typename vtype2::reg_t,
    // NumPy✔️❌:           typename opmask_t = typename vtype1::opmask_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE reg_t1 cmp_merge(reg_t1 in1,
    // NumPy✔️❌:                                       reg_t1 in2,
    // NumPy✔️❌:                                       reg_t2 &indexes1,
    // NumPy✔️❌:                                       reg_t2 indexes2,
    // NumPy✔️❌:                                       opmask_t mask)
    // NumPy✔️❌: {
    // NumPy✔️❌:     reg_t1 tmp_keys = cmp_merge<vtype1>(in1, in2, mask);
    // NumPy✔️❌:     indexes1 = vtype2::mask_mov(
    // NumPy✔️❌:             indexes2,
    // NumPy✔️❌:             resize_mask<vtype1, vtype2>(vtype1::eq(tmp_keys, in1)),
    // NumPy✔️❌:             indexes1);
    // NumPy✔️❌:     return tmp_keys; // 0 -> min, 1 -> max
    // NumPy✔️❌: }
    let first = *values;
    for i in 0..8 {
        let key = if mask & (1 << i) == 0 {
            if other[i].key < first[i].key {
                other[i].key
            } else {
                first[i].key
            }
        } else if other[i].key > first[i].key {
            other[i].key
        } else {
            first[i].key
        };
        values[i] = Lane {
            key,
            index: if key == first[i].key {
                first[i].index
            } else {
                other[i].index
            },
        };
    }
}
fn permutation(values: Register, width: usize, reverse: bool) -> Register {
    // reverse_n reverses within each width; swap_n swaps its two halves.
    // The masks are compile-time constants in the source register networks.
    let xor = if reverse { width - 1 } else { width / 2 };
    std::array::from_fn(|i| values[i ^ xor])
}
fn sort_register(values: &mut Register) {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename vtype1,
    // NumPy✔️❌:           typename vtype2,
    // NumPy✔️❌:           typename reg_t = typename vtype1::reg_t,
    // NumPy✔️❌:           typename index_type = typename vtype2::reg_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE reg_t sort_reg_8lanes(reg_t key_reg, index_type &index_reg)
    // NumPy✔️❌: {
    // NumPy✔️❌:     using key_swizzle = typename vtype1::swizzle_ops;
    // NumPy✔️❌:     using index_swizzle = typename vtype2::swizzle_ops;
    // NumPy✔️❌:
    // NumPy✔️❌:     const auto oxAA = convert_int_to_mask<vtype1>(0xAA);
    // NumPy✔️❌:     const auto oxCC = convert_int_to_mask<vtype1>(0xCC);
    // NumPy✔️❌:     const auto oxF0 = convert_int_to_mask<vtype1>(0xF0);
    // NumPy✔️❌:
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template reverse_n<vtype1, 2>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template reverse_n<vtype2, 2>(index_reg),
    // NumPy✔️❌:             oxAA);
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template reverse_n<vtype1, 4>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template reverse_n<vtype2, 4>(index_reg),
    // NumPy✔️❌:             oxCC);
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template swap_n<vtype1, 2>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template swap_n<vtype2, 2>(index_reg),
    // NumPy✔️❌:             oxAA);
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template reverse_n<vtype1, 8>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template reverse_n<vtype2, 8>(index_reg),
    // NumPy✔️❌:             oxF0);
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template swap_n<vtype1, 4>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template swap_n<vtype2, 4>(index_reg),
    // NumPy✔️❌:             oxCC);
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template swap_n<vtype1, 2>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template swap_n<vtype2, 2>(index_reg),
    // NumPy✔️❌:             oxAA);
    // NumPy✔️❌:     return key_reg;
    // NumPy✔️❌: }
    for (width, reverse, mask) in [
        (2, true, 0xAA),
        (4, true, 0xCC),
        (2, false, 0xAA),
        (8, true, 0xF0),
        (4, false, 0xCC),
        (2, false, 0xAA),
    ] {
        let other = permutation(*values, width, reverse);
        compare_merge(values, other, mask);
    }
}
fn merge_register(values: &mut Register) {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename vtype1,
    // NumPy✔️❌:           typename vtype2,
    // NumPy✔️❌:           typename reg_t = typename vtype1::reg_t,
    // NumPy✔️❌:           typename index_type = typename vtype2::reg_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE reg_t bitonic_merge_reg_8lanes(reg_t key_reg,
    // NumPy✔️❌:                                                     index_type &index_reg)
    // NumPy✔️❌: {
    // NumPy✔️❌:     using key_swizzle = typename vtype1::swizzle_ops;
    // NumPy✔️❌:     using index_swizzle = typename vtype2::swizzle_ops;
    // NumPy✔️❌:
    // NumPy✔️❌:     const auto oxAA = convert_int_to_mask<vtype1>(0xAA);
    // NumPy✔️❌:     const auto oxCC = convert_int_to_mask<vtype1>(0xCC);
    // NumPy✔️❌:     const auto oxF0 = convert_int_to_mask<vtype1>(0xF0);
    // NumPy✔️❌:
    // NumPy✔️❌:     // 1) half_cleaner[8]: compare 0-4, 1-5, 2-6, 3-7
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template swap_n<vtype1, 8>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template swap_n<vtype2, 8>(index_reg),
    // NumPy✔️❌:             oxF0);
    // NumPy✔️❌:     // 2) half_cleaner[4]
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template swap_n<vtype1, 4>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template swap_n<vtype2, 4>(index_reg),
    // NumPy✔️❌:             oxCC);
    // NumPy✔️❌:     // 3) half_cleaner[1]
    // NumPy✔️❌:     key_reg = cmp_merge<vtype1, vtype2>(
    // NumPy✔️❌:             key_reg,
    // NumPy✔️❌:             key_swizzle::template swap_n<vtype1, 2>(key_reg),
    // NumPy✔️❌:             index_reg,
    // NumPy✔️❌:             index_swizzle::template swap_n<vtype2, 2>(index_reg),
    // NumPy✔️❌:             oxAA);
    // NumPy✔️❌:     return key_reg;
    // NumPy✔️❌: }
    for (width, mask) in [(8, 0xF0), (4, 0xCC), (2, 0xAA)] {
        let other = permutation(*values, width, false);
        compare_merge(values, other, mask);
    }
}
fn compare_exchange(values: &mut [Register], left: usize, right: usize) {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename vtype1,
    // NumPy✔️❌:           typename vtype2,
    // NumPy✔️❌:           typename reg_t1 = typename vtype1::reg_t,
    // NumPy✔️❌:           typename reg_t2 = typename vtype2::reg_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE void
    // NumPy✔️❌: COEX(reg_t1 &key1, reg_t1 &key2, reg_t2 &index1, reg_t2 &index2)
    // NumPy✔️❌: {
    // NumPy✔️❌:     reg_t1 key_t1 = vtype1::min(key1, key2);
    // NumPy✔️❌:     reg_t1 key_t2 = vtype1::max(key1, key2);
    // NumPy✔️❌:
    // NumPy✔️❌:     auto eqMask = resize_mask<vtype1, vtype2>(vtype1::eq(key_t1, key1));
    // NumPy✔️❌:
    // NumPy✔️❌:     reg_t2 index_t1 = vtype2::mask_mov(index2, eqMask, index1);
    // NumPy✔️❌:     reg_t2 index_t2 = vtype2::mask_mov(index1, eqMask, index2);
    // NumPy✔️❌:
    // NumPy✔️❌:     key1 = key_t1;
    // NumPy✔️❌:     key2 = key_t2;
    // NumPy✔️❌:     index1 = index_t1;
    // NumPy✔️❌:     index2 = index_t2;
    // NumPy✔️❌: }
    let a = values[left];
    let b = values[right];
    for i in 0..8 {
        if b[i].key < a[i].key {
            values[left][i] = b[i];
            values[right][i] = a[i];
        }
    }
}
fn merge_vectors(values: &mut [Register]) {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename keyType, typename valueType, int numVecs>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE void bitonic_merge_n_vec(typename keyType::reg_t *keys,
    // NumPy✔️❌:                                               typename valueType::reg_t *values)
    // NumPy✔️❌: {
    // NumPy✔️❌:     // Do the reverse part
    // NumPy✔️❌:     if constexpr (numVecs == 2) {
    // NumPy✔️❌:         keys[1] = keyType::reverse(keys[1]);
    // NumPy✔️❌:         values[1] = valueType::reverse(values[1]);
    // NumPy✔️❌:         COEX<keyType, valueType>(keys[0], keys[1], values[0], values[1]);
    // NumPy✔️❌:         keys[1] = keyType::reverse(keys[1]);
    // NumPy✔️❌:         values[1] = valueType::reverse(values[1]);
    // NumPy✔️❌:     }
    // NumPy✔️❌:     else if constexpr (numVecs > 2) {
    // NumPy✔️❌:         // Reverse upper half
    // NumPy✔️❌:         X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:         for (int i = 0; i < numVecs / 2; i++) {
    // NumPy✔️❌:             keys[numVecs - i - 1] = keyType::reverse(keys[numVecs - i - 1]);
    // NumPy✔️❌:             values[numVecs - i - 1]
    // NumPy✔️❌:                     = valueType::reverse(values[numVecs - i - 1]);
    // NumPy✔️❌:
    // NumPy✔️❌:             COEX<keyType, valueType>(keys[i],
    // NumPy✔️❌:                                      keys[numVecs - i - 1],
    // NumPy✔️❌:                                      values[i],
    // NumPy✔️❌:                                      values[numVecs - i - 1]);
    // NumPy✔️❌:
    // NumPy✔️❌:             keys[numVecs - i - 1] = keyType::reverse(keys[numVecs - i - 1]);
    // NumPy✔️❌:             values[numVecs - i - 1]
    // NumPy✔️❌:                     = valueType::reverse(values[numVecs - i - 1]);
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌:     // Call cleaner
    // NumPy✔️❌:     bitonic_clean_n_vec<keyType, valueType, numVecs>(keys, values);
    // NumPy✔️❌:
    // NumPy✔️❌:     // Now do bitonic_merge
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:     for (int i = 0; i < numVecs; i++) {
    // NumPy✔️❌:         bitonic_merge_dispatch<keyType, valueType>(keys[i], values[i]);
    // NumPy✔️❌:     }
    // NumPy✔️❌: }
    let count = values.len();
    for i in 0..count / 2 {
        let right = count - i - 1;
        values[right].reverse();
        compare_exchange(values, i, right);
        values[right].reverse();
    }
    let mut width = count / 2;
    while width >= 2 {
        for j in (0..count).step_by(width) {
            for i in 0..width / 2 {
                compare_exchange(values, j + i, j + i + width / 2);
            }
        }
        width /= 2;
    }
    for value in values {
        merge_register(value);
    }
}
fn network(keys: &[f64], indices: &mut [usize]) {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename keyType, typename indexType, int numVecs>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE void
    // NumPy✔️❌: argsort_n_vec(const typename keyType::type_t *keys, arrsize_t *indices, int N)
    // NumPy✔️❌: {
    // NumPy✔️❌:     using kreg_t = typename keyType::reg_t;
    // NumPy✔️❌:     using ireg_t = typename indexType::reg_t;
    // NumPy✔️❌:
    // NumPy✔️❌:     static_assert(numVecs > 0, "numVecs should be > 0");
    // NumPy✔️❌:     if constexpr (numVecs > 1) {
    // NumPy✔️❌:         if (N * 2 <= numVecs * keyType::numlanes) {
    // NumPy✔️❌:             argsort_n_vec<keyType, indexType, numVecs / 2>(keys, indices, N);
    // NumPy✔️❌:             return;
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌:     kreg_t keyVecs[numVecs];
    // NumPy✔️❌:     ireg_t indexVecs[numVecs];
    // NumPy✔️❌:
    // NumPy✔️❌:     // Generate masks for loading and storing
    // NumPy✔️❌:     typename keyType::opmask_t ioMasks[numVecs - numVecs / 2];
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:     for (int i = numVecs / 2, j = 0; i < numVecs; i++, j++) {
    // NumPy✔️❌:         uint64_t num_to_read
    // NumPy✔️❌:                 = std::min((uint64_t)std::max(0, N - i * keyType::numlanes),
    // NumPy✔️❌:                            (uint64_t)keyType::numlanes);
    // NumPy✔️❌:         ioMasks[j] = keyType::get_partial_loadmask(num_to_read);
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌:     // Unmasked part of the load
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:     for (int i = 0; i < numVecs / 2; i++) {
    // NumPy✔️❌:         indexVecs[i] = indexType::loadu(indices + i * indexType::numlanes);
    // NumPy✔️❌:         keyVecs[i]
    // NumPy✔️❌:                 = keyType::i64gather(keys, indices + i * indexType::numlanes);
    // NumPy✔️❌:     }
    // NumPy✔️❌:     // Masked part of the load
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:     for (int i = numVecs / 2; i < numVecs; i++) {
    // NumPy✔️❌:         indexVecs[i] = indexType::mask_loadu(
    // NumPy✔️❌:                 indexType::zmm_max(),
    // NumPy✔️❌:                 resize_mask<keyType, indexType>(ioMasks[i - numVecs / 2]),
    // NumPy✔️❌:                 indices + i * indexType::numlanes);
    // NumPy✔️❌:
    // NumPy✔️❌:         keyVecs[i] = keyType::template mask_i64gather<sizeof(
    // NumPy✔️❌:                 typename keyType::type_t)>(keyType::zmm_max(),
    // NumPy✔️❌:                                            ioMasks[i - numVecs / 2],
    // NumPy✔️❌:                                            indexVecs[i],
    // NumPy✔️❌:                                            keys);
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌:     // Sort each loaded vector
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:     for (int i = 0; i < numVecs; i++) {
    // NumPy✔️❌:         sort_vec_dispatch<keyType, indexType>(keyVecs[i], indexVecs[i]);
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌:     // Run the full merger
    // NumPy✔️❌:     bitonic_fullmerge_n_vec<keyType, indexType, numVecs>(keyVecs, indexVecs);
    // NumPy✔️❌:
    // NumPy✔️❌:     // Unmasked part of the store
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:     for (int i = 0; i < numVecs / 2; i++) {
    // NumPy✔️❌:         indexType::storeu(indices + i * indexType::numlanes, indexVecs[i]);
    // NumPy✔️❌:     }
    // NumPy✔️❌:     // Masked part of the store
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(64)
    // NumPy✔️❌:     for (int i = numVecs / 2, j = 0; i < numVecs; i++, j++) {
    // NumPy✔️❌:         indexType::mask_storeu(
    // NumPy✔️❌:                 indices + i * indexType::numlanes,
    // NumPy✔️❌:                 resize_mask<keyType, indexType>(ioMasks[i - numVecs / 2]),
    // NumPy✔️❌:                 indexVecs[i]);
    // NumPy✔️❌:     }
    // NumPy✔️❌: }
    let count = indices.len().div_ceil(8).next_power_of_two();
    let mut values = [[EMPTY; 8]; 32];
    for (i, &index) in indices.iter().enumerate() {
        values[i / 8][i % 8] = Lane {
            key: keys[index],
            index,
        };
    }
    for value in &mut values[..count] {
        sort_register(value);
    }
    let mut width = 2;
    while width <= count {
        for chunk in values[..count].chunks_exact_mut(width) {
            merge_vectors(chunk);
        }
        width *= 2;
    }
    for (i, index) in indices.iter_mut().enumerate() {
        *index = values[i / 8][i % 8].index;
    }
}
fn store_partition(
    keys: &[f64],
    indices: &mut [usize],
    left: usize,
    right: usize,
    values: [usize; 8],
    pivot: f64,
    smallest: &mut f64,
    biggest: &mut f64,
) -> usize {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename vtype,
    // NumPy✔️❌:           typename argtype,
    // NumPy✔️❌:           typename type_t,
    // NumPy✔️❌:           typename reg_t,
    // NumPy✔️❌:           typename argreg_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE int32_t partition_vec_avx512(type_t *arg,
    // NumPy✔️❌:                                                   arrsize_t left,
    // NumPy✔️❌:                                                   arrsize_t right,
    // NumPy✔️❌:                                                   const argreg_t arg_vec,
    // NumPy✔️❌:                                                   const reg_t curr_vec,
    // NumPy✔️❌:                                                   const reg_t pivot_vec,
    // NumPy✔️❌:                                                   reg_t *smallest_vec,
    // NumPy✔️❌:                                                   reg_t *biggest_vec)
    // NumPy✔️❌: {
    // NumPy✔️❌:     /* which elements are larger than the pivot */
    // NumPy✔️❌:     typename vtype::opmask_t gt_mask = vtype::ge(curr_vec, pivot_vec);
    // NumPy✔️❌:     int32_t amount_gt_pivot = _mm_popcnt_u32((int32_t)gt_mask);
    // NumPy✔️❌:     argtype::mask_compressstoreu(
    // NumPy✔️❌:             arg + left, vtype::knot_opmask(gt_mask), arg_vec);
    // NumPy✔️❌:     argtype::mask_compressstoreu(
    // NumPy✔️❌:             arg + right - amount_gt_pivot, gt_mask, arg_vec);
    // NumPy✔️❌:     *smallest_vec = vtype::min(curr_vec, *smallest_vec);
    // NumPy✔️❌:     *biggest_vec = vtype::max(curr_vec, *biggest_vec);
    // NumPy✔️❌:     return amount_gt_pivot;
    // NumPy✔️❌: }
    let greater = values.iter().filter(|&&i| keys[i] >= pivot).count();
    let mut lo = left;
    let mut hi = right - greater;
    for index in values {
        let value = keys[index];
        if value < *smallest {
            *smallest = value;
        }
        if value > *biggest {
            *biggest = value;
        }
        if value >= pivot {
            indices[hi] = index;
            hi += 1;
        } else {
            indices[lo] = index;
            lo += 1;
        }
    }
    greater
}
fn partition(
    keys: &[f64],
    indices: &mut [usize],
    mut left: usize,
    mut right: usize,
    pivot: f64,
) -> (usize, f64, f64) {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename vtype,
    // NumPy✔️❌:           typename argtype,
    // NumPy✔️❌:           int num_unroll,
    // NumPy✔️❌:           typename type_t = typename vtype::type_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE arrsize_t argpartition_unrolled(const type_t *arr,
    // NumPy✔️❌:                                                      arrsize_t *arg,
    // NumPy✔️❌:                                                      arrsize_t left,
    // NumPy✔️❌:                                                      arrsize_t right,
    // NumPy✔️❌:                                                      type_t pivot,
    // NumPy✔️❌:                                                      type_t *smallest,
    // NumPy✔️❌:                                                      type_t *biggest)
    // NumPy✔️❌: {
    // NumPy✔️❌:     if (right - left <= 8 * num_unroll * vtype::numlanes) {
    // NumPy✔️❌:         return argpartition<vtype, argtype>(
    // NumPy✔️❌:                 arr, arg, left, right, pivot, smallest, biggest);
    // NumPy✔️❌:     }
    // NumPy✔️❌:     /* make array length divisible by vtype::numlanes , shortening the array */
    // NumPy✔️❌:     for (int32_t i = ((right - left) % (num_unroll * vtype::numlanes)); i > 0;
    // NumPy✔️❌:          --i) {
    // NumPy✔️❌:         *smallest = std::min(*smallest, arr[arg[left]], comparison_func<vtype>);
    // NumPy✔️❌:         *biggest = std::max(*biggest, arr[arg[left]], comparison_func<vtype>);
    // NumPy✔️❌:         if (!comparison_func<vtype>(arr[arg[left]], pivot)) {
    // NumPy✔️❌:             std::swap(arg[left], arg[--right]);
    // NumPy✔️❌:         }
    // NumPy✔️❌:         else {
    // NumPy✔️❌:             ++left;
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌:     if (left == right)
    // NumPy✔️❌:         return left; /* less than vtype::numlanes elements in the array */
    // NumPy✔️❌:
    // NumPy✔️❌:     using reg_t = typename vtype::reg_t;
    // NumPy✔️❌:     using argreg_t = typename argtype::reg_t;
    // NumPy✔️❌:     reg_t pivot_vec = vtype::set1(pivot);
    // NumPy✔️❌:     reg_t min_vec = vtype::set1(*smallest);
    // NumPy✔️❌:     reg_t max_vec = vtype::set1(*biggest);
    // NumPy✔️❌:
    // NumPy✔️❌:     // first and last vtype::numlanes values are partitioned at the end
    // NumPy✔️❌:     reg_t vec_left[num_unroll], vec_right[num_unroll];
    // NumPy✔️❌:     argreg_t argvec_left[num_unroll], argvec_right[num_unroll];
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(8)
    // NumPy✔️❌:     for (int ii = 0; ii < num_unroll; ++ii) {
    // NumPy✔️❌:         argvec_left[ii] = argtype::loadu(arg + left + vtype::numlanes * ii);
    // NumPy✔️❌:         vec_left[ii] = vtype::i64gather(arr, arg + left + vtype::numlanes * ii);
    // NumPy✔️❌:         argvec_right[ii] = argtype::loadu(
    // NumPy✔️❌:                 arg + (right - vtype::numlanes * (num_unroll - ii)));
    // NumPy✔️❌:         vec_right[ii] = vtype::i64gather(
    // NumPy✔️❌:                 arr, arg + (right - vtype::numlanes * (num_unroll - ii)));
    // NumPy✔️❌:     }
    // NumPy✔️❌:     // store points of the vectors
    // NumPy✔️❌:     arrsize_t r_store = right - vtype::numlanes;
    // NumPy✔️❌:     arrsize_t l_store = left;
    // NumPy✔️❌:     // indices for loading the elements
    // NumPy✔️❌:     left += num_unroll * vtype::numlanes;
    // NumPy✔️❌:     right -= num_unroll * vtype::numlanes;
    // NumPy✔️❌:     while (right - left != 0) {
    // NumPy✔️❌:         argreg_t arg_vec[num_unroll];
    // NumPy✔️❌:         reg_t curr_vec[num_unroll];
    // NumPy✔️❌:         /*
    // NumPy✔️❌:          * if fewer elements are stored on the right side of the array,
    // NumPy✔️❌:          * then next elements are loaded from the right side,
    // NumPy✔️❌:          * otherwise from the left side
    // NumPy✔️❌:          */
    // NumPy✔️❌:         if ((r_store + vtype::numlanes) - right < left - l_store) {
    // NumPy✔️❌:             right -= num_unroll * vtype::numlanes;
    // NumPy✔️❌:             X86_SIMD_SORT_UNROLL_LOOP(8)
    // NumPy✔️❌:             for (int ii = 0; ii < num_unroll; ++ii) {
    // NumPy✔️❌:                 arg_vec[ii]
    // NumPy✔️❌:                         = argtype::loadu(arg + right + ii * vtype::numlanes);
    // NumPy✔️❌:                 curr_vec[ii] = vtype::i64gather(
    // NumPy✔️❌:                         arr, arg + right + ii * vtype::numlanes);
    // NumPy✔️❌:             }
    // NumPy✔️❌:         }
    // NumPy✔️❌:         else {
    // NumPy✔️❌:             X86_SIMD_SORT_UNROLL_LOOP(8)
    // NumPy✔️❌:             for (int ii = 0; ii < num_unroll; ++ii) {
    // NumPy✔️❌:                 arg_vec[ii] = argtype::loadu(arg + left + ii * vtype::numlanes);
    // NumPy✔️❌:                 curr_vec[ii] = vtype::i64gather(
    // NumPy✔️❌:                         arr, arg + left + ii * vtype::numlanes);
    // NumPy✔️❌:             }
    // NumPy✔️❌:             left += num_unroll * vtype::numlanes;
    // NumPy✔️❌:         }
    // NumPy✔️❌:         // partition the current vector and save it on both sides of the array
    // NumPy✔️❌:         X86_SIMD_SORT_UNROLL_LOOP(8)
    // NumPy✔️❌:         for (int ii = 0; ii < num_unroll; ++ii) {
    // NumPy✔️❌:             int32_t amount_gt_pivot
    // NumPy✔️❌:                     = partition_vec<vtype, argtype>(arg,
    // NumPy✔️❌:                                                     l_store,
    // NumPy✔️❌:                                                     r_store + vtype::numlanes,
    // NumPy✔️❌:                                                     arg_vec[ii],
    // NumPy✔️❌:                                                     curr_vec[ii],
    // NumPy✔️❌:                                                     pivot_vec,
    // NumPy✔️❌:                                                     &min_vec,
    // NumPy✔️❌:                                                     &max_vec);
    // NumPy✔️❌:             l_store += (vtype::numlanes - amount_gt_pivot);
    // NumPy✔️❌:             r_store -= amount_gt_pivot;
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌:     /* partition and save vec_left and vec_right */
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(8)
    // NumPy✔️❌:     for (int ii = 0; ii < num_unroll; ++ii) {
    // NumPy✔️❌:         int32_t amount_gt_pivot
    // NumPy✔️❌:                 = partition_vec<vtype, argtype>(arg,
    // NumPy✔️❌:                                                 l_store,
    // NumPy✔️❌:                                                 r_store + vtype::numlanes,
    // NumPy✔️❌:                                                 argvec_left[ii],
    // NumPy✔️❌:                                                 vec_left[ii],
    // NumPy✔️❌:                                                 pivot_vec,
    // NumPy✔️❌:                                                 &min_vec,
    // NumPy✔️❌:                                                 &max_vec);
    // NumPy✔️❌:         l_store += (vtype::numlanes - amount_gt_pivot);
    // NumPy✔️❌:         r_store -= amount_gt_pivot;
    // NumPy✔️❌:     }
    // NumPy✔️❌:     X86_SIMD_SORT_UNROLL_LOOP(8)
    // NumPy✔️❌:     for (int ii = 0; ii < num_unroll; ++ii) {
    // NumPy✔️❌:         int32_t amount_gt_pivot
    // NumPy✔️❌:                 = partition_vec<vtype, argtype>(arg,
    // NumPy✔️❌:                                                 l_store,
    // NumPy✔️❌:                                                 r_store + vtype::numlanes,
    // NumPy✔️❌:                                                 argvec_right[ii],
    // NumPy✔️❌:                                                 vec_right[ii],
    // NumPy✔️❌:                                                 pivot_vec,
    // NumPy✔️❌:                                                 &min_vec,
    // NumPy✔️❌:                                                 &max_vec);
    // NumPy✔️❌:         l_store += (vtype::numlanes - amount_gt_pivot);
    // NumPy✔️❌:         r_store -= amount_gt_pivot;
    // NumPy✔️❌:     }
    // NumPy✔️❌:     *smallest = vtype::reducemin(min_vec);
    // NumPy✔️❌:     *biggest = vtype::reducemax(max_vec);
    // NumPy✔️❌:     return l_store;
    // NumPy✔️❌: }
    // argsort_ calls this fixed unroll=4 path only for >256 elements;
    // its source <=256 argpartition branch is unreachable after network dispatch.
    let mut smallest = f64::INFINITY;
    let mut biggest = f64::NEG_INFINITY;
    for _ in 0..(right - left) % 32 {
        let value = keys[indices[left]];
        if value < smallest {
            smallest = value;
        }
        if value > biggest {
            biggest = value;
        }
        if value >= pivot {
            right -= 1;
            indices.swap(left, right);
        } else {
            left += 1;
        }
    }
    if left == right {
        return (left, smallest, biggest);
    }
    let first: [[usize; 8]; 4] = std::array::from_fn(|i| {
        indices[left + i * 8..left + (i + 1) * 8]
            .try_into()
            .unwrap()
    });
    let last: [[usize; 8]; 4] = std::array::from_fn(|i| {
        indices[right - 32 + i * 8..right - 24 + i * 8]
            .try_into()
            .unwrap()
    });
    let mut l_store = left;
    let mut r_store = right as isize - 8;
    left += 32;
    right -= 32;
    while right != left {
        let values: [[usize; 8]; 4];
        if (r_store + 8) - (right as isize) < (left as isize) - (l_store as isize) {
            right -= 32;
            values = std::array::from_fn(|i| {
                indices[right + i * 8..right + (i + 1) * 8]
                    .try_into()
                    .unwrap()
            });
        } else {
            values = std::array::from_fn(|i| {
                indices[left + i * 8..left + (i + 1) * 8]
                    .try_into()
                    .unwrap()
            });
            left += 32;
        }
        for value in values {
            let greater = store_partition(
                keys,
                indices,
                l_store,
                (r_store + 8) as usize,
                value,
                pivot,
                &mut smallest,
                &mut biggest,
            );
            l_store += 8 - greater;
            r_store -= greater as isize;
        }
    }
    for values in [first, last] {
        for value in values {
            let greater = store_partition(
                keys,
                indices,
                l_store,
                (r_store + 8) as usize,
                value,
                pivot,
                &mut smallest,
                &mut biggest,
            );
            l_store += 8 - greater;
            r_store -= greater as isize;
        }
    }
    (l_store, smallest, biggest)
}
fn pivot(keys: &[f64], indices: &[usize]) -> f64 {
    // Scalar lane emulation keeps source order and allocation/complexity bounds,
    // but costs more scalar operations/register copies than source AVX512.

    // NumPy✔️❌: template <typename vtype, typename type_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE type_t get_pivot_64bit(const type_t *arr,
    // NumPy✔️❌:                                             arrsize_t *arg,
    // NumPy✔️❌:                                             const arrsize_t left,
    // NumPy✔️❌:                                             const arrsize_t right)
    // NumPy✔️❌: {
    // NumPy✔️❌:     if constexpr (vtype::numlanes == 8) {
    // NumPy✔️❌:         if (right - left >= vtype::numlanes) {
    // NumPy✔️❌:             // median of 8
    // NumPy✔️❌:             arrsize_t size = (right - left) / 8;
    // NumPy✔️❌:             using reg_t = typename vtype::reg_t;
    // NumPy✔️❌:             reg_t rand_vec = vtype::set(arr[arg[left + size]],
    // NumPy✔️❌:                                         arr[arg[left + 2 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 3 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 4 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 5 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 6 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 7 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 8 * size]]);
    // NumPy✔️❌:             // pivot will never be a nan, since there are no nan's!
    // NumPy✔️❌:             reg_t sort = vtype::sort_vec(rand_vec);
    // NumPy✔️❌:             return ((type_t *)&sort)[4];
    // NumPy✔️❌:         }
    // NumPy✔️❌:         else {
    // NumPy✔️❌:             return arr[arg[right]];
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌:     else if constexpr (vtype::numlanes == 4) {
    // NumPy✔️❌:         if (right - left >= vtype::numlanes) {
    // NumPy✔️❌:             // median of 4
    // NumPy✔️❌:             arrsize_t size = (right - left) / 4;
    // NumPy✔️❌:             using reg_t = typename vtype::reg_t;
    // NumPy✔️❌:             reg_t rand_vec = vtype::set(arr[arg[left + size]],
    // NumPy✔️❌:                                         arr[arg[left + 2 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 3 * size]],
    // NumPy✔️❌:                                         arr[arg[left + 4 * size]]);
    // NumPy✔️❌:             // pivot will never be a nan, since there are no nan's!
    // NumPy✔️❌:             reg_t sort = vtype::sort_vec(rand_vec);
    // NumPy✔️❌:             return ((type_t *)&sort)[2];
    // NumPy✔️❌:         }
    // NumPy✔️❌:         else {
    // NumPy✔️❌:             return arr[arg[right]];
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌: }
    let step = (indices.len() - 1) / 8;
    // vtype::set lists high lane first. Pivot key, not its index, is observed.
    let mut values = std::array::from_fn(|i| Lane {
        key: keys[indices[(8 - i) * step]],
        index: 0,
    });
    sort_register(&mut values);
    values[4].key
}
fn recurse(keys: &[f64], indices: &mut [usize], depth: usize) {
    // NumPy✔️❌: template <typename vtype, typename argtype, typename type_t>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE void argsort_(const type_t *arr,
    // NumPy✔️❌:                                    arrsize_t *arg,
    // NumPy✔️❌:                                    arrsize_t left,
    // NumPy✔️❌:                                    arrsize_t right,
    // NumPy✔️❌:                                    arrsize_t max_iters,
    // NumPy✔️❌:                                    arrsize_t task_threshold)
    // NumPy✔️❌: {
    // NumPy✔️❌:     /*
    // NumPy✔️❌:      * Resort to std::sort if quicksort isnt making any progress
    // NumPy✔️❌:      */
    // NumPy✔️✔️:     if (max_iters <= 0) {
    // NumPy✔️✔️:         std_argsort(arr, arg, left, right + 1);
    // NumPy✔️✔️:         return;
    // NumPy✔️✔️:     }
    // NumPy✔️❌:     /*
    // NumPy✔️❌:      * Base case: use bitonic networks to sort arrays <= 64
    // NumPy✔️❌:      */
    // NumPy✔️❌:     if (right + 1 - left <= 256) {
    // NumPy✔️❌:         argsort_n<vtype, argtype, 256>(
    // NumPy✔️❌:                 arr, arg + left, (int32_t)(right + 1 - left));
    // NumPy✔️❌:         return;
    // NumPy✔️❌:     }
    // NumPy✔️❌:     type_t pivot = get_pivot_64bit<vtype>(arr, arg, left, right);
    // NumPy✔️❌:     type_t smallest = vtype::type_max();
    // NumPy✔️❌:     type_t biggest = vtype::type_min();
    // NumPy✔️❌:     arrsize_t pivot_index = argpartition_unrolled<vtype, argtype, 4>(
    // NumPy✔️❌:             arr, arg, left, right + 1, pivot, &smallest, &biggest);
    // NumPy✔️❌: #ifdef XSS_COMPILE_OPENMP
    // NumPy✔️❌:     if (pivot != smallest) {
    // NumPy✔️❌:         bool parallel_left = (pivot_index - left) > task_threshold;
    // NumPy✔️❌:         if (parallel_left) {
    // NumPy✔️❌: #pragma omp task
    // NumPy✔️❌:             argsort_<vtype, argtype>(arr,
    // NumPy✔️❌:                                      arg,
    // NumPy✔️❌:                                      left,
    // NumPy✔️❌:                                      pivot_index - 1,
    // NumPy✔️❌:                                      max_iters - 1,
    // NumPy✔️❌:                                      task_threshold);
    // NumPy✔️❌:         }
    // NumPy✔️❌:         else {
    // NumPy✔️❌:             argsort_<vtype, argtype>(arr,
    // NumPy✔️❌:                                      arg,
    // NumPy✔️❌:                                      left,
    // NumPy✔️❌:                                      pivot_index - 1,
    // NumPy✔️❌:                                      max_iters - 1,
    // NumPy✔️❌:                                      task_threshold);
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌:     if (pivot != biggest) {
    // NumPy✔️❌:         bool parallel_right = (right - pivot_index) > task_threshold;
    // NumPy✔️❌:
    // NumPy✔️❌:         if (parallel_right) {
    // NumPy✔️❌: #pragma omp task
    // NumPy✔️❌:             argsort_<vtype, argtype>(arr,
    // NumPy✔️❌:                                      arg,
    // NumPy✔️❌:                                      pivot_index,
    // NumPy✔️❌:                                      right,
    // NumPy✔️❌:                                      max_iters - 1,
    // NumPy✔️❌:                                      task_threshold);
    // NumPy✔️❌:         }
    // NumPy✔️❌:         else {
    // NumPy✔️❌:             argsort_<vtype, argtype>(arr,
    // NumPy✔️❌:                                      arg,
    // NumPy✔️❌:                                      pivot_index,
    // NumPy✔️❌:                                      right,
    // NumPy✔️❌:                                      max_iters - 1,
    // NumPy✔️❌:                                      task_threshold);
    // NumPy✔️❌:         }
    // NumPy✔️❌:     }
    // NumPy✔️❌: #else
    // NumPy✔️❌:     UNUSED(task_threshold);
    // NumPy✔️❌:     if (pivot != smallest)
    // NumPy✔️❌:         argsort_<vtype, argtype>(
    // NumPy✔️❌:                 arr, arg, left, pivot_index - 1, max_iters - 1, 0);
    // NumPy✔️❌:     if (pivot != biggest)
    // NumPy✔️❌:         argsort_<vtype, argtype>(
    // NumPy✔️❌:                 arr, arg, pivot_index, right, max_iters - 1, 0);
    // NumPy✔️❌: #endif
    // NumPy✔️❌: }
    if depth == 0 {
        cosmolkit_core::source_sort_by(indices, |a, b| keys[*a] < keys[*b]);
        return;
    }
    if indices.len() <= 256 {
        network(keys, indices);
        return;
    }
    let value = pivot(keys, indices);
    let (cut, smallest, biggest) = partition(keys, indices, 0, indices.len(), value);
    let (left, right) = indices.split_at_mut(cut);
    if value != smallest {
        recurse(keys, left, depth - 1);
    }
    if value != biggest {
        recurse(keys, right, depth - 1);
    }
}
/// The named NumPy reference's ascending float64 argsort; keys are untouched.
pub(super) fn argsort(keys: &[f64]) -> Vec<usize> {
    // NumPy✔️❌: template <typename...>
    // NumPy✔️❌:           typename half_vector>
    // NumPy✔️❌: X86_SIMD_SORT_INLINE void xss_argsort(const T *arr,
    // NumPy✔️❌:                                       arrsize_t *arg,
    // NumPy✔️❌:                                       arrsize_t arrsize,
    // NumPy✔️❌:                                       bool hasnan = false,
    // NumPy✔️❌:                                       bool descending = false,
    // NumPy✔️❌:                                       bool nans_last = true)
    // NumPy✔️❌: {
    // NumPy✔️❌:
    // NumPy✔️❌:     using vectype = typename std::conditional<sizeof(T) == sizeof(int32_t),
    // NumPy✔️❌:                                               half_vector<T>,
    // NumPy✔️❌:                                               full_vector<T>>::type;
    // NumPy✔️❌:
    // NumPy✔️❌:     using argtype =
    // NumPy✔️❌:             typename std::conditional<sizeof(arrsize_t) == sizeof(int32_t),
    // NumPy✔️❌:                                       half_vector<arrsize_t>,
    // NumPy✔️❌:                                       full_vector<arrsize_t>>::type;
    // NumPy✔️❌:
    // NumPy✔️❌:     if (arrsize > 1) {
    // NumPy✔️❌:         /* simdargsort does not work for float/double arrays with nan */
    // NumPy✔️❌:         if constexpr (xss::fp::is_floating_point_v<T>) {
    // NumPy✔️✔️:             if ((hasnan) && (array_has_nan<vectype>(arr, arrsize))) {
    // NumPy✔️✔️:                 std_argsort_withnan(
    // NumPy✔️✔️:                         arr, arg, 0, arrsize, nans_last, descending);
    // NumPy✔️✔️:                 return;
    // NumPy✔️✔️:             }
    // NumPy✔️❌:         }
    // NumPy✔️❌:         UNUSED(hasnan);
    // NumPy✔️❌:         UNUSED(nans_last);
    // NumPy✔️❌:
    // NumPy✔️❌:         /* early exit for already sorted arrays: float/double with nan never reach here*/
    // NumPy✔️❌:         auto comp = descending ? Comparator<vectype, true>::STDSortComparator
    // NumPy✔️❌:                                : Comparator<vectype, false>::STDSortComparator;
    // NumPy✔️✔️:         if (std::is_sorted(arr, arr + arrsize, comp)) { return; }
    // NumPy✔️❌:
    // NumPy✔️❌: #ifdef XSS_COMPILE_OPENMP
    // NumPy✔️❌:
    // NumPy✔️❌:         bool use_parallel = arrsize > 10000;
    // NumPy✔️❌:
    // NumPy✔️❌:         if (use_parallel) {
    // NumPy✔️❌:             int thread_count = xss_get_num_threads();
    // NumPy✔️❌:             arrsize_t task_threshold
    // NumPy✔️❌:                     = std::max((arrsize_t)10000, arrsize / 100);
    // NumPy✔️❌:
    // NumPy✔️❌:             // We use omp parallel and then omp single to setup the threads that will run the omp task calls in qsort_
    // NumPy✔️❌:             // The omp single prevents multiple threads from running the initial qsort_ simultaneously and causing problems
    // NumPy✔️❌:             // Note that we do not use the if(...) clause built into OpenMP, because it causes a performance regression for small arrays
    // NumPy✔️❌: #pragma omp parallel num_threads(thread_count)
    // NumPy✔️❌: #pragma omp single
    // NumPy✔️❌:             argsort_<vectype, argtype>(arr,
    // NumPy✔️❌:                                        arg,
    // NumPy✔️❌:                                        0,
    // NumPy✔️❌:                                        arrsize - 1,
    // NumPy✔️❌:                                        2 * (arrsize_t)log2(arrsize),
    // NumPy✔️❌:                                        task_threshold);
    // NumPy✔️❌: #pragma omp taskwait
    // NumPy✔️❌:         }
    // NumPy✔️❌:         else {
    // NumPy✔️❌:             argsort_<vectype, argtype>(arr,
    // NumPy✔️❌:                                        arg,
    // NumPy✔️❌:                                        0,
    // NumPy✔️❌:                                        arrsize - 1,
    // NumPy✔️❌:                                        2 * (arrsize_t)log2(arrsize),
    // NumPy✔️❌:                                        std::numeric_limits<arrsize_t>::max());
    // NumPy✔️❌:         }
    // NumPy✔️❌: #else
    // NumPy✔️❌:         argsort_<vectype, argtype>(
    // NumPy✔️❌:                 arr, arg, 0, arrsize - 1, 2 * (arrsize_t)log2(arrsize), 0);
    // NumPy✔️❌: #endif
    // NumPy✔️❌:
    // NumPy✔️❌:         if (descending) { std::reverse(arg, arg + arrsize); }
    // NumPy✔️❌:     }
    // NumPy✔️❌:
    // NumPy✔️❌: #ifdef __MMX__
    // NumPy✔️❌:     // Workaround for compiler bug generating MMX instructions without emms
    // NumPy✔️❌:     _mm_empty();
    // NumPy✔️❌: #endif
    // NumPy✔️❌: }
    let mut indices: Vec<_> = (0..keys.len()).collect();
    if keys.len() > 1 {
        if keys.iter().any(|v| v.is_nan()) {
            cosmolkit_core::source_sort_by(&mut indices, |a, b| {
                let (a, b) = (keys[*a], keys[*b]);
                if a.is_nan() {
                    false
                } else if b.is_nan() {
                    true
                } else {
                    a < b
                }
            });
        } else {
            // Source std::is_sorted returns before any partitioning. Apart
            // from its O(n) fast path, this preserves equal-key indices in an
            // already sorted input beyond the 256-element network boundary.
            if keys.windows(2).all(|pair| !(pair[1] < pair[0])) {
                return indices;
            }
            let depth = 2 * ((keys.len() as f64).log2() as usize);
            recurse(keys, &mut indices, depth);
        }
    }
    indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numpy_avx512_argsort_keeps_source_ties_and_partition_order() {
        // Fixed NumPy 2.5.3/X86_V4 reference results, including unstable ties.
        assert_eq!(
            argsort(&[1., 0., 1., 0., 1., 0., 1., 0.]),
            [1, 3, 7, 5, 2, 0, 4, 6]
        );
        let keys: Vec<_> = (0..17).map(|i| ((i * 7) % 5) as f64).collect();
        assert_eq!(
            argsort(&keys),
            [0, 5, 15, 10, 8, 13, 3, 1, 11, 6, 16, 4, 14, 9, 7, 2, 12]
        );
        assert_eq!(argsort(&[f64::INFINITY, 1.0, f64::INFINITY]), [1, 0, 2]);
        assert_eq!(argsort(&[f64::NAN, 1.0, 0.0, f64::NAN]), [2, 1, 0, 3]);
        // The first input beyond the source's network threshold must use the
        // same complete partition path. Unique integer keys give an exact,
        // independently determined expected permutation without another sorter.
        let keys: Vec<_> = (0..257).map(|i| ((i * 149) % 257) as f64).collect();
        let expected: Vec<_> = (0..257)
            .map(|key| keys.iter().position(|&v| v == key as f64).unwrap())
            .collect();
        assert_eq!(argsort(&keys), expected);
        // xss_argsort's source early exit applies before its 256-element
        // network/partition dispatch, including already ordered equal keys.
        let equal = vec![1.0; 257];
        assert_eq!(argsort(&equal), (0..equal.len()).collect::<Vec<_>>());
        let ordered: Vec<_> = (0..257).map(|i| (i / 5) as f64).collect();
        assert_eq!(argsort(&ordered), (0..ordered.len()).collect::<Vec<_>>());
        assert!(argsort(&[]).is_empty());
        assert_eq!(argsort(&[0.]), [0]);
    }
}
