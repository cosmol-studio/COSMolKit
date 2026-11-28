//! GraphDescriptors.py/Graphs.py information content, RDKit 2026.09.1.
//! Copyright (C) 2001-2022 Greg Landrum and other RDKit contributors (BSD).
//! NumPy reduction and OpenBLAS operation-order portions retain their notices in LICENSE.
use cosmolkit_model::TopologyBlock;
use std::fmt;
mod entropy_log;

/// A supplied square matrix is borrowed; every entry equal to one is an adjacency entry.
#[derive(Debug, Clone, Copy)]
pub enum IpcMatrixInput<'a> {
    /// Checked square row storage.
    Rows(&'a [Vec<f64>]),
    /// A dense square matrix in its native storage.
    Dense(&'a cosmolkit_core::DenseMatrix),
}

/// Input shape failure for an Ipc distance-matrix override.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum IpcError {
    /// A matrix row count or row width did not match the molecule's atom count.
    InvalidDistanceMatrixShape {
        /// Row with the wrong width, or `None` for the overall row count/dense dimension.
        row: Option<usize>,
        /// Observed row count, width, or dense dimension.
        actual: usize,
        /// Required atom count.
        expected: usize,
    },
    /// A default topology bond referred to an atom outside the stored atom rows.
    InvalidBondEndpoint {
        /// Bond row.
        bond: usize,
        /// Endpoint name.
        endpoint: &'static str,
        /// Stored atom index.
        atom: usize,
        /// Number of stored atoms.
        atom_count: usize,
    },
}

impl fmt::Display for IpcError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDistanceMatrixShape {
                row: Some(row),
                actual,
                expected,
            } => write!(
                f,
                "Ipc distance matrix row {row} has {actual} columns; expected {expected}"
            ),
            Self::InvalidDistanceMatrixShape {
                row: None,
                actual,
                expected,
            } => write!(
                f,
                "Ipc distance matrix has dimension {actual}; expected {expected}"
            ),
            Self::InvalidBondEndpoint {
                bond,
                endpoint,
                atom,
                atom_count,
            } => write!(
                f,
                "Ipc topology bond {bond} {endpoint} endpoint {atom} is outside {atom_count} atoms"
            ),
        }
    }
}

impl std::error::Error for IpcError {}

/// Computes total information content from the supplied graph.
pub fn ipc(topology: &TopologyBlock) -> Result<f64, IpcError> {
    let graph = graph_from_topology(topology)?;
    Ok(information_content(&graph, false))
}

/// Computes total information content using a checked, borrowed matrix override.
pub fn ipc_with_distance_matrix(
    topology: &TopologyBlock,
    matrix: IpcMatrixInput<'_>,
) -> Result<f64, IpcError> {
    let graph = graph_from_matrix(topology.atoms.len(), matrix)?;
    Ok(information_content(&graph, false))
}

/// Computes average information content from the supplied graph.
pub fn avg_ipc(topology: &TopologyBlock) -> Result<f64, IpcError> {
    let graph = graph_from_topology(topology)?;
    Ok(information_content(&graph, true))
}

/// Computes average information content using a checked, borrowed matrix override.
pub fn avg_ipc_with_distance_matrix(
    topology: &TopologyBlock,
    matrix: IpcMatrixInput<'_>,
) -> Result<f64, IpcError> {
    let graph = graph_from_matrix(topology.atoms.len(), matrix)?;
    Ok(information_content(&graph, true))
}

struct BinaryGraph {
    dimension: usize,
    adjacency: Vec<f64>,
}

fn graph_from_topology(topology: &TopologyBlock) -> Result<BinaryGraph, IpcError> {
    // RDKit✔️🔝: dMat = Chem.GetDistanceMatrix(mol, False)
    // RDKit✔️🔝: adjMat = numpy.equal(dMat, 1)
    // Matrices.cpp, the useBO=false path:
    // RDKit✔️🔝:   } else {
    // RDKit✔️🔝:     contrib = 1.0;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   dMat[i * nAts + j] = contrib;
    // RDKit✔️🔝:   dMat[j * nAts + i] = contrib;
    // Unit-length edges are exactly the entries at distance one: no
    // Floyd-Warshall work is needed to recover these entries. Duplicate
    // bonds assign the same entry instead of accumulating their weights.
    let dimension = topology.atoms.len();
    let mut adjacency = vec![0.0; dimension * dimension];
    for (bond_index, bond) in topology.bonds.iter().enumerate() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        for (endpoint, atom) in [("begin", begin), ("end", end)] {
            if atom >= dimension {
                return Err(IpcError::InvalidBondEndpoint {
                    bond: bond_index,
                    endpoint,
                    atom,
                    atom_count: dimension,
                });
            }
        }
        adjacency[begin * dimension + end] = 1.0;
        adjacency[end * dimension + begin] = 1.0;
    }
    Ok(BinaryGraph {
        dimension,
        adjacency,
    })
}

fn graph_from_matrix(
    dimension: usize,
    matrix: IpcMatrixInput<'_>,
) -> Result<BinaryGraph, IpcError> {
    match matrix {
        IpcMatrixInput::Rows(rows) => graph_from_rows(dimension, rows.iter().map(Vec::as_slice)),
        IpcMatrixInput::Dense(matrix) => {
            graph_from_rows(dimension, matrix.values().chunks(matrix.dimension().max(1)))
        }
    }
}

fn graph_from_rows<'a>(
    dimension: usize,
    rows: impl ExactSizeIterator<Item = &'a [f64]>,
) -> Result<BinaryGraph, IpcError> {
    if rows.len() != dimension {
        return Err(IpcError::InvalidDistanceMatrixShape {
            row: None,
            actual: rows.len(),
            expected: dimension,
        });
    }
    let mut adjacency = Vec::with_capacity(dimension * dimension);
    for (row, values) in rows.enumerate() {
        if values.len() != dimension {
            return Err(IpcError::InvalidDistanceMatrixShape {
                row: Some(row),
                actual: values.len(),
                expected: dimension,
            });
        }
        // RDKit✔️✔️: adjMat = numpy.equal(dMat, 1)
        adjacency.extend(values.iter().map(|&value| f64::from(value == 1.0)));
    }
    Ok(BinaryGraph {
        dimension,
        adjacency,
    })
}

fn information_content(graph: &BinaryGraph, average: bool) -> f64 {
    // RDKit❗✔️: adjMat = numpy.equal(dMat, 1)
    // RDKit❗✔️: cPoly = abs(Graphs.CharacteristicPolynomial(mol, adjMat))
    // RDKit❗✔️: if avg:
    // RDKit❗✔️:   return entropy.InfoEntropy(cPoly)
    // RDKit❗✔️: else:
    // RDKit❗✔️:   return sum(cPoly) * entropy.InfoEntropy(cPoly)
    let mut coefficients = characteristic_polynomial(graph);
    for value in &mut coefficients {
        *value = value.abs();
    }
    let entropy = info_entropy(&coefficients);
    if average {
        entropy
    } else {
        coefficients.iter().sum::<f64>() * entropy
    }
}

fn characteristic_polynomial(graph: &BinaryGraph) -> Vec<f64> {
    // RDKit❗✔️: nAtoms = mol.GetNumAtoms()
    // RDKit❗✔️: A = mat
    // RDKit❗✔️: I = 1. * numpy.identity(nAtoms)
    // RDKit❗✔️: An = A
    // RDKit❗✔️: res = numpy.zeros(nAtoms + 1, float)
    // RDKit❗✔️: res[0] = 1.0
    // RDKit❗✔️: for n in range(1, nAtoms + 1):
    // RDKit❗✔️:   res[n] = 1. / n * numpy.trace(An)
    // RDKit❗✔️:   Bn = An - res[n] * I
    // RDKit❗✔️:   An = numpy.dot(A, Bn)
    // RDKit❗✔️: res[1:] *= -1
    // RDKit❗✔️: return res
    let n = graph.dimension;
    let mut result = vec![0.0; n + 1];
    result[0] = 1.0;
    if n == 0 {
        return result;
    }
    let mut power = graph.adjacency.clone();
    let mut next = vec![0.0; n * n];
    // Reuse the two n² buffers; no per-iteration identity/matrix allocation.
    for coefficient in result.iter_mut().skip(1).enumerate() {
        let (index, value) = coefficient;
        *value = (1.0 / (index + 1) as f64) * trace_sum(&power, n, 0, n);
        for (i, row) in power.chunks_exact_mut(n).enumerate() {
            for (j, entry) in row.iter_mut().enumerate() {
                // Multiplication by an identity off-diagonal must remain:
                // Inf*0 or NaN*0 is NaN in NumPy, not a no-op.
                *entry -= *value * f64::from(i == j);
            }
        }
        matrix_product(&graph.adjacency, &power, &mut next, n);
        std::mem::swap(&mut power, &mut next);
    }
    for value in &mut result[1..] {
        *value *= -1.0;
    }
    result
}

fn trace_sum(matrix: &[f64], stride: usize, begin: usize, len: usize) -> f64 {
    // NumPy 2.5.3 _core/src/umath/loops_utils.h.src (PW_BLOCKSIZE=128):
    // NumPy✔️✔️: static inline @type@
    // NumPy✔️✔️: @TYPE@_pairwise_sum(char *a, npy_intp n, npy_intp stride)
    // NumPy✔️✔️: {
    // NumPy✔️✔️:     if (n < 8) {
    // NumPy✔️✔️:         npy_intp i;
    // NumPy✔️✔️:         /*
    // NumPy✔️✔️:          * Start with -0 to preserve -0 values.  The reason is that summing
    // NumPy✔️✔️:          * only -0 should return -0, but `0 + -0 == 0` while `-0 + -0 == -0`.
    // NumPy✔️✔️:          */
    // NumPy✔️✔️:         @type@ res = -0.0;
    // NumPy✔️✔️:
    // NumPy✔️✔️:         for (i = 0; i < n; i++) {
    // NumPy✔️✔️:             res += @trf@(*((@dtype@*)(a + i * stride)));
    // NumPy✔️✔️:         }
    // NumPy✔️✔️:         return res;
    // NumPy✔️✔️:     }
    // NumPy✔️✔️:     else if (n <= PW_BLOCKSIZE) {
    // NumPy✔️✔️:         npy_intp i;
    // NumPy✔️✔️:         @type@ r[8], res;
    // NumPy✔️✔️:
    // NumPy✔️✔️:         /*
    // NumPy✔️✔️:          * sum a block with 8 accumulators
    // NumPy✔️✔️:          * 8 times unroll reduces blocksize to 16 and allows vectorization with
    // NumPy✔️✔️:          * avx without changing summation ordering
    // NumPy✔️✔️:          */
    // NumPy✔️✔️:         r[0] = @trf@(*((@dtype@ *)(a + 0 * stride)));
    // NumPy✔️✔️:         r[1] = @trf@(*((@dtype@ *)(a + 1 * stride)));
    // NumPy✔️✔️:         r[2] = @trf@(*((@dtype@ *)(a + 2 * stride)));
    // NumPy✔️✔️:         r[3] = @trf@(*((@dtype@ *)(a + 3 * stride)));
    // NumPy✔️✔️:         r[4] = @trf@(*((@dtype@ *)(a + 4 * stride)));
    // NumPy✔️✔️:         r[5] = @trf@(*((@dtype@ *)(a + 5 * stride)));
    // NumPy✔️✔️:         r[6] = @trf@(*((@dtype@ *)(a + 6 * stride)));
    // NumPy✔️✔️:         r[7] = @trf@(*((@dtype@ *)(a + 7 * stride)));
    // NumPy✔️✔️:
    // NumPy✔️✔️:         for (i = 8; i < n - (n % 8); i += 8) {
    // NumPy✔️✔️:             /* small blocksizes seems to mess with hardware prefetch */
    // NumPy✔️✔️:             NPY_PREFETCH(a + (i + 512/(npy_intp)sizeof(@dtype@))*stride, 0, 3);
    // NumPy✔️✔️:             r[0] += @trf@(*((@dtype@ *)(a + (i + 0) * stride)));
    // NumPy✔️✔️:             r[1] += @trf@(*((@dtype@ *)(a + (i + 1) * stride)));
    // NumPy✔️✔️:             r[2] += @trf@(*((@dtype@ *)(a + (i + 2) * stride)));
    // NumPy✔️✔️:             r[3] += @trf@(*((@dtype@ *)(a + (i + 3) * stride)));
    // NumPy✔️✔️:             r[4] += @trf@(*((@dtype@ *)(a + (i + 4) * stride)));
    // NumPy✔️✔️:             r[5] += @trf@(*((@dtype@ *)(a + (i + 5) * stride)));
    // NumPy✔️✔️:             r[6] += @trf@(*((@dtype@ *)(a + (i + 6) * stride)));
    // NumPy✔️✔️:             r[7] += @trf@(*((@dtype@ *)(a + (i + 7) * stride)));
    // NumPy✔️✔️:         }
    // NumPy✔️✔️:
    // NumPy✔️✔️:         /* accumulate now to avoid stack spills for single peel loop */
    // NumPy✔️✔️:         res = ((r[0] + r[1]) + (r[2] + r[3])) +
    // NumPy✔️✔️:               ((r[4] + r[5]) + (r[6] + r[7]));
    // NumPy✔️✔️:
    // NumPy✔️✔️:         /* do non multiple of 8 rest */
    // NumPy✔️✔️:         for (; i < n; i++) {
    // NumPy✔️✔️:             res += @trf@(*((@dtype@ *)(a + i * stride)));
    // NumPy✔️✔️:         }
    // NumPy✔️✔️:         return res;
    // NumPy✔️✔️:     }
    // NumPy✔️✔️:     else {
    // NumPy✔️✔️:         /* divide by two but avoid non-multiples of unroll factor */
    // NumPy✔️✔️:         npy_intp n2 = n / 2;
    // NumPy✔️✔️:
    // NumPy✔️✔️:         n2 -= n2 % 8;
    // NumPy✔️✔️:         return @TYPE@_pairwise_sum(a, n2, stride) +
    // NumPy✔️✔️:                @TYPE@_pairwise_sum(a + n2 * stride, n - n2, stride);
    // NumPy✔️✔️:     }
    // NumPy✔️✔️: }
    let at = |i| matrix[(begin + i) * (stride + 1)];
    if len < 8 {
        return (0..len).fold(-0.0, |sum, i| sum + at(i));
    }
    if len > 128 {
        let half = (len / 2) / 8 * 8;
        return trace_sum(matrix, stride, begin, half)
            + trace_sum(matrix, stride, begin + half, len - half);
    }
    let mut lanes = std::array::from_fn::<_, 8, _>(at);
    let end = len / 8 * 8;
    for i in (8..end).step_by(8) {
        for lane in 0..8 {
            lanes[lane] += at(i + lane);
        }
    }
    let mut sum = ((lanes[0] + lanes[1]) + (lanes[2] + lanes[3]))
        + ((lanes[4] + lanes[5]) + (lanes[6] + lanes[7]));
    for i in end..len {
        sum += at(i);
    }
    // add.reduce starts with its +0 identity before the pairwise block.
    0.0 + sum
}

fn matrix_product(left: &[f64], right: &[f64], output: &mut [f64], n: usize) {
    // OpenBLAS v0.3.34 kernel/x86_64/dgemm_small_kernel_permit_skylakex.c:
    // OpenBLAS❗❗: if (MNK > 100.0*100.0*100.0)  // disable for big size matrix
    // OpenBLAS❗❗:     return 0;
    // dgemm_small_kernel_nn_skylakex.c (row-major CBLAS swaps A and B):
    // OpenBLAS❗❗: for (; i < m8; i += 8) {
    // OpenBLAS❗❗: int mm = M - i;
    // OpenBLAS❗❗: if (!mm) return 0;
    // OpenBLAS❗❗: if (mm > 4 || K < 16) {
    // OpenBLAS❗❗:     for (k = 0; k < K; k++) {
    // OpenBLAS❗❗: } else {
    // OpenBLAS❗❗:     BLASLONG k8 = K & ~7;
    // OpenBLAS❗❗:     for (k = 0; k < k8; k += 8) {
    // OpenBLAS❗❗: #define MATMUL_512(M, N) result##M##N = _mm512_fmadd_pd(Aval##M, Bval##N, result##M##N)
    // One factor is binary adjacency (0 or 1), so its multiplication of a
    // finite value is exact. Separate multiply/add has the same rounding as
    // FMA and avoids a software FMA call for every matrix element on targets
    // without hardware FMA. Keep every zero product for nonfinite propagation.
    // Fixed reference arithmetic is independent of the executing CPU's SIMD width.
    // Bulk outputs accumulate k in order; the small-kernel tail uses eight lanes.
    if n > 100 {
        blocked_matrix_product(left, right, output, n);
        return;
    }
    let remainder = n % 8;
    let packet_tail = (16..=100).contains(&n) && (1..=4).contains(&remainder);
    let bulk = if packet_tail { n - remainder } else { n };
    output.fill(0.0);
    for (i, row) in output.chunks_exact_mut(n).enumerate() {
        for k in 0..n {
            let factor = left[i * n + k];
            for (entry, &value) in row[..bulk].iter_mut().zip(&right[k * n..k * n + bulk]) {
                *entry += factor * value;
            }
        }
        for (j, entry) in row.iter_mut().enumerate().skip(bulk) {
            let mut lanes = [0.0; 8];
            for k in 0..n {
                lanes[k % 8] += left[i * n + k] * right[k * n + j];
            }
            *entry = if remainder == 4 || i < n / 4 * 4 {
                ((lanes[0] + lanes[1]) + (lanes[2] + lanes[3]))
                    + ((lanes[4] + lanes[5]) + (lanes[6] + lanes[7]))
            } else {
                ((lanes[6] + lanes[2]) + (lanes[4] + lanes[0]))
                    + ((lanes[7] + lanes[3]) + (lanes[5] + lanes[1]))
            };
        }
    }
}

fn blocked_matrix_product(left: &[f64], right: &[f64], output: &mut [f64], n: usize) {
    // OpenBLAS v0.3.34 param.h, SKYLAKEX (the pinned NumPy wheel):
    // OpenBLAS✔️❗: #define DGEMM_DEFAULT_P 192
    // OpenBLAS✔️❗: #define DGEMM_DEFAULT_Q 384
    // OpenBLAS✔️❗: #define DGEMM_DEFAULT_R 8640
    // driver/level3/level3.c:
    // OpenBLAS✔️❗:       min_i = m_to - m_from;
    // OpenBLAS✔️❗:       if (min_i >= GEMM_P * 2) {
    // OpenBLAS✔️❗: 	min_i = GEMM_P;
    // OpenBLAS✔️❗:       } else {
    // OpenBLAS✔️❗: 	if (min_i > GEMM_P) {
    // OpenBLAS✔️❗: 	  min_i = ((min_i / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
    // The row-major BLAS call swaps operands: source m is our column axis.
    // Its initial m block is called with short n panels; subsequent m blocks
    // use the full n panel and its 12-column kernels. All m-block boundaries
    // are multiples of 16, so only the final n%8 columns can use K lanes.
    let first_columns = if n >= 384 {
        192
    } else if n > 192 {
        (n / 2).div_ceil(16) * 16
    } else {
        n
    };
    let tail_start = n - n % 8;
    let two_lane_end = tail_start + if n % 8 >= 4 { 4 } else { 0 };
    // R=8640 is divisible by 12, so panel-local full n12 blocks end here.
    let full_rows = n / 12 * 12;
    let mut sums = vec![0.0; n];
    output.fill(0.0);
    let mut begin = 0;
    while begin < n {
        // OpenBLAS✔️❗:       min_l = k - ls;
        // OpenBLAS✔️❗:       if (min_l >= GEMM_Q * 2) {
        // OpenBLAS✔️❗: 	min_l  = GEMM_Q;
        // OpenBLAS✔️❗:       } else {
        // OpenBLAS✔️❗: 	if (min_l > GEMM_Q) {
        // OpenBLAS✔️❗: 	  min_l = ((min_l / 2 + GEMM_UNROLL_M - 1)/GEMM_UNROLL_M) * GEMM_UNROLL_M;
        let remaining = n - begin;
        let width = if remaining >= 768 {
            384
        } else if remaining > 384 {
            (remaining / 2).div_ceil(16) * 16
        } else {
            remaining
        };
        let end = begin + width;
        for (i, row) in output.chunks_exact_mut(n).enumerate() {
            let split = if i < full_rows && tail_start >= first_columns {
                tail_start
            } else {
                n
            };
            sums[..split].fill(0.0);
            for k in begin..end {
                let factor = left[i * n + k];
                for (sum, &value) in sums[..split].iter_mut().zip(&right[k * n..k * n + split]) {
                    *sum += factor * value;
                }
            }
            for (entry, &sum) in row.iter_mut().zip(&sums[..split]) {
                *entry += sum;
            }
            for (j, entry) in row.iter_mut().enumerate().skip(split) {
                // kernel/x86_64/dgemm_kernel_16x2_skylakex.c:
                // KERNEL_k1m4n12 consumes pairs, then unit_merge_to_ymm:
                // OpenBLAS✔️❗:   "vextractf64x4 $1,%%zmm"#c1_no",%%ymm30; vaddpd %%ymm"#c1_no",%%ymm30,%%ymm"#c1_no";"
                // KERNEL_k1m2n12 / KERNEL_k1m1n12 consume four K entries,
                // then unit_merge_to_xmm combines adjacent lane pairs:
                // OpenBLAS✔️❗:   "vaddpd %%xmm20,%%xmm21,%%xmm20; vaddpd %%xmm22,%%xmm23,%%xmm22; vaddpd %%xmm20,%%xmm22,%%xmm"#c1_no";"
                // Remaining K entries are added after this merge. Source
                // SAVE adds each completed K block to C, not to its lanes.
                let count = if j < two_lane_end { 2 } else { 4 };
                let stop = begin + width / count * count;
                let mut lanes = [0.0; 4];
                for k in begin..stop {
                    lanes[(k - begin) % count] += left[i * n + k] * right[k * n + j];
                }
                let mut sum = if count == 2 {
                    lanes[0] + lanes[1]
                } else {
                    (lanes[0] + lanes[1]) + (lanes[2] + lanes[3])
                };
                for k in stop..end {
                    sum += left[i * n + k] * right[k * n + j];
                }
                *entry += sum;
            }
        }
        begin = end;
    }
}

fn info_entropy(values: &[f64]) -> f64 {
    // RDKit❗✔️: template <class T>
    // RDKit❗✔️: double InfoEntropy(const T *tPtr, long int dim) {
    // RDKit❗✔️:   int i;
    // RDKit❗✔️:   T nInstances = 0;
    // RDKit❗✔️:   double accum = 0.0, d;
    // RDKit❗✔️:   for (i = 0; i < dim; i++) {
    // RDKit❗✔️:     nInstances += tPtr[i];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (nInstances != 0) {
    // RDKit❗✔️:     for (i = 0; i < dim; i++) {
    // RDKit❗✔️:       d = (double)tPtr[i] / nInstances;
    // RDKit❗✔️:       if (d != 0) {
    // RDKit❗✔️:         accum += -d * log(d);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return accum / log(2.0);
    // RDKit❗✔️: }
    let mut total = 0.0;
    for &value in values {
        total += value;
    }
    if total == 0.0 {
        return 0.0;
    }

    let mut accum = 0.0;
    for &value in values {
        let probability = value / total;
        if probability != 0.0 {
            accum += -probability * entropy_log::probability_ln(probability);
        }
    }
    accum / std::f64::consts::LN_2
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn large_products_preserve_kernel_lanes_and_k_blocks() {
        // OpenBLAS 0.3.34 / NumPy 2.5.3: the same four products reduce to
        // 1 sequentially, 2 in two lanes and 0 in four lanes. Exercise all
        // m4/m2/m1 tails and the final non-n12 row without a molecular fixture.
        for n in 193..200 {
            let mut left = vec![0.0; n * n];
            let mut right = vec![0.0; n * n];
            let mut output = vec![0.0; n * n];
            for i in 0..n {
                left[i * n..i * n + 4].fill(1.0);
            }
            for (k, value) in [2_f64.powi(54), 1.0, -2_f64.powi(54), 1.0]
                .into_iter()
                .enumerate()
            {
                right[k * n..(k + 1) * n].fill(value);
            }
            matrix_product(&left, &right, &mut output, n);
            let start = n - n % 8;
            for j in 0..n {
                let expected: f64 = if j < start {
                    1.0
                } else if n % 8 >= 4 && j < start + 4 {
                    2.0
                } else {
                    0.0
                };
                assert_eq!(output[j].to_bits(), expected.to_bits(), "n={n}, j={j}");
                // The final source n block has fewer than twelve columns.
                assert_eq!(output[(n - 1) * n + j].to_bits(), 1.0_f64.to_bits());
            }
        }
        let n = 385;
        let mut left = vec![0.0; n * n];
        let mut right = vec![0.0; n * n];
        for (k, value) in [
            (0, 2_f64.powi(54)),
            (1, 1.0),
            (192, -2_f64.powi(54)),
            (193, 1.0),
        ] {
            left[k] = 1.0;
            right[k * n] = value;
        }
        let mut output = vec![0.0; n * n];
        matrix_product(&left, &right, &mut output, n);
        assert_eq!(output[0].to_bits(), 0.0_f64.to_bits());
    }
    fn graph(n: usize, edges: &[(usize, usize)]) -> BinaryGraph {
        let mut adjacency = vec![0.0; n * n];
        for &(i, j) in edges {
            adjacency[i * n + j] = 1.0;
            adjacency[j * n + i] = 1.0;
        }
        BinaryGraph {
            dimension: n,
            adjacency,
        }
    }
    #[test]
    fn coefficients_follow_the_released_float_recurrence() {
        let triangle = characteristic_polynomial(&graph(3, &[(0, 1), (1, 2), (0, 2)]));
        assert_eq!(triangle, [1.0, -0.0, -3.0, -2.0]);
        assert_eq!(
            characteristic_polynomial(&graph(4, &[(0, 1), (2, 3)])),
            [1.0, -0.0, -2.0, -0.0, 1.0]
        );
        assert_eq!(information_content(&graph(0, &[]), false), 0.0);
        assert_eq!(information_content(&graph(4, &[]), true), 0.0);
    }
    #[test]
    fn asymmetric_overrides_use_the_whole_matrix() {
        let directed = graph_from_rows(
            2,
            [vec![0.0, 0.0], vec![1.0, 0.0]].iter().map(Vec::as_slice),
        )
        .unwrap();
        assert_eq!(information_content(&directed, false), 0.0);
        let looped = graph_from_rows(1, [vec![1.0]].iter().map(Vec::as_slice)).unwrap();
        assert_eq!(information_content(&looped, false), 2.0);
        assert_eq!(information_content(&looped, true), 1.0);
    }
}
