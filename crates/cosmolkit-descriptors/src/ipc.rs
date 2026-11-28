//! GraphDescriptors.py/Graphs.py information content, RDKit 2026.09.2.
//! Copyright (C) 2001-2022 Greg Landrum and other RDKit contributors (BSD).
//! NumPy reduction and OpenBLAS operation-order portions retain their notices in LICENSE.
use cosmolkit_model::TopologyBlock;
use std::fmt;
mod entropy_log;
mod exact;
mod lapack_blas;
mod lapack_blas_level3;
mod lapack_general;
mod lapack_symmetric;
mod lapack_tridiagonal;
mod numpy_argsort;
mod numpy_convolve;

/// Characteristic-polynomial algorithm used for molecular information content.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum IpcMethod {
    /// Preserve the pinned RDKit 2026.09.2 eigenvalue expansion.
    #[default]
    RdkitCompatible,
    /// RDKit PR #9884 exact integer coefficients followed by binary64 entropy.
    ExactPolynomial,
}

/// Total information content using an explicitly selected polynomial algorithm.
pub fn ipc_with_method(topology: &TopologyBlock, method: IpcMethod) -> Result<f64, IpcError> {
    match method {
        IpcMethod::RdkitCompatible => ipc(topology),
        IpcMethod::ExactPolynomial => exact::calculate(topology, false),
    }
}

/// Average information content using an explicitly selected polynomial algorithm.
pub fn avg_ipc_with_method(topology: &TopologyBlock, method: IpcMethod) -> Result<f64, IpcError> {
    match method {
        IpcMethod::RdkitCompatible => avg_ipc(topology),
        IpcMethod::ExactPolynomial => exact::calculate(topology, true),
    }
}

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
    /// LAPACK returned positive INFO instead of a complete eigenvalue result.
    EigenvaluesDidNotConverge {
        /// The source routine reporting incomplete convergence.
        routine: &'static str,
        /// Original positive LAPACK INFO, with that routine's meaning.
        info: usize,
    },
    /// A numerical routine rejected its argument at the source boundary.
    InvalidNumericalArgument {
        /// Original source routine.
        routine: &'static str,
        /// One-based source argument position.
        argument: u8,
    },
    /// A checked numerical view exceeded its supplied backing storage.
    InsufficientNumericalStorage {
        /// Original source routine.
        routine: &'static str,
        /// Source matrix or vector name.
        matrix: &'static str,
        /// Minimum required number of backing elements.
        required: usize,
        /// Number of supplied backing elements.
        actual: usize,
    },
    /// Numerical dimensions exceeded the machine index range.
    NumericalDimensionOverflow,
    /// The source DLASRT fixed sorting stack could not hold another interval.
    NumericalSortingStackExhausted,
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
            Self::EigenvaluesDidNotConverge { routine, info } => write!(
                f,
                "Ipc eigenvalues did not converge in {routine} (INFO={info})"
            ),
            Self::InvalidNumericalArgument { routine, argument } => {
                write!(f, "Ipc {routine} argument {argument} is invalid")
            }
            Self::InsufficientNumericalStorage {
                routine,
                matrix,
                required,
                actual,
            } => write!(
                f,
                "Ipc {routine} {matrix} requires {required} elements, got {actual}"
            ),
            Self::NumericalDimensionOverflow => {
                f.write_str("Ipc numerical dimensions overflow the index range")
            }
            Self::NumericalSortingStackExhausted => {
                f.write_str("Ipc DLASRT fixed sorting stack exhausted")
            }
        }
    }
}

impl std::error::Error for IpcError {}

impl From<lapack_blas_level3::Blas3Error> for IpcError {
    fn from(error: lapack_blas_level3::Blas3Error) -> Self {
        use lapack_blas_level3::Blas3Error as E;
        match error {
            E::InvalidArgument { routine, argument } => {
                Self::InvalidNumericalArgument { routine, argument }
            }
            E::StorageTooShort {
                routine,
                matrix,
                required,
                actual,
            } => Self::InsufficientNumericalStorage {
                routine,
                matrix,
                required,
                actual,
            },
            E::DimensionOverflow => Self::NumericalDimensionOverflow,
        }
    }
}

impl From<lapack_tridiagonal::DsterfError> for IpcError {
    fn from(error: lapack_tridiagonal::DsterfError) -> Self {
        use lapack_tridiagonal::DsterfError as E;
        match error {
            E::InvalidOffDiagonalLength { expected, actual } => {
                Self::InsufficientNumericalStorage {
                    routine: "DSTERF",
                    matrix: "E",
                    required: expected,
                    actual,
                }
            }
            E::DimensionOverflow => Self::NumericalDimensionOverflow,
            E::InvalidArgument { routine, argument } => {
                Self::InvalidNumericalArgument { routine, argument }
            }
            E::SortingStackExhausted => Self::NumericalSortingStackExhausted,
            E::NonConvergence { remaining } => Self::EigenvaluesDidNotConverge {
                routine: "DSTERF",
                info: remaining,
            },
        }
    }
}

impl From<lapack_symmetric::SymmetricError> for IpcError {
    fn from(error: lapack_symmetric::SymmetricError) -> Self {
        match error {
            lapack_symmetric::SymmetricError::Tridiagonal(error) => error.into(),
            lapack_symmetric::SymmetricError::Blas(error) => error.into(),
        }
    }
}

impl From<lapack_general::GeneralError> for IpcError {
    fn from(error: lapack_general::GeneralError) -> Self {
        match error {
            lapack_general::GeneralError::InvalidArgument { routine, argument } => {
                Self::InvalidNumericalArgument { routine, argument }
            }
            lapack_general::GeneralError::NonConvergence { info } => {
                Self::EigenvaluesDidNotConverge {
                    routine: "DHSEQR",
                    info,
                }
            }
            lapack_general::GeneralError::Blas(error) => error.into(),
            lapack_general::GeneralError::Scalar(error) => error.into(),
            lapack_general::GeneralError::DimensionOverflow => Self::NumericalDimensionOverflow,
        }
    }
}

/// Computes total information content from the supplied graph.
pub fn ipc(topology: &TopologyBlock) -> Result<f64, IpcError> {
    let graph = graph_from_topology(topology)?;
    information_content(&graph, false)
}

/// Computes total information content using a checked, borrowed matrix override.
pub fn ipc_with_distance_matrix(
    topology: &TopologyBlock,
    matrix: IpcMatrixInput<'_>,
) -> Result<f64, IpcError> {
    let graph = graph_from_matrix(topology.atoms.len(), matrix)?;
    information_content(&graph, false)
}

/// Computes average information content from the supplied graph.
pub fn avg_ipc(topology: &TopologyBlock) -> Result<f64, IpcError> {
    let graph = graph_from_topology(topology)?;
    information_content(&graph, true)
}

/// Computes average information content using a checked, borrowed matrix override.
pub fn avg_ipc_with_distance_matrix(
    topology: &TopologyBlock,
    matrix: IpcMatrixInput<'_>,
) -> Result<f64, IpcError> {
    let graph = graph_from_matrix(topology.atoms.len(), matrix)?;
    information_content(&graph, true)
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

fn information_content(graph: &BinaryGraph, average: bool) -> Result<f64, IpcError> {
    // RDKit❗✔️: adjMat = numpy.equal(dMat, 1)
    // RDKit❗✔️: cPoly = abs(Graphs.CharacteristicPolynomial(mol, adjMat))
    // RDKit❗✔️: if avg:
    // RDKit❗✔️:   return entropy.InfoEntropy(cPoly)
    // RDKit❗✔️: else:
    // RDKit❗✔️:   return sum(cPoly) * entropy.InfoEntropy(cPoly)
    let mut coefficients = characteristic_polynomial(graph)?;
    for value in &mut coefficients {
        *value = value.abs();
    }
    let entropy = info_entropy(&coefficients);
    Ok(if average {
        entropy
    } else {
        coefficients.iter().sum::<f64>() * entropy
    })
}

fn characteristic_polynomial(graph: &BinaryGraph) -> Result<Vec<f64>, IpcError> {
    // RDKit❗✔️:   A = numpy.asarray(A, dtype=float)
    // RDKit❗✔️:   # A molecular adjacency matrix is symmetric, and eigvalsh is backward stable for that case.
    // RDKit❗✔️:   # The test must be EXACT: eigvalsh reads only one triangle, so a merely near-symmetric matrix
    // RDKit❗✔️:   # would silently get the polynomial of its symmetrised version instead of its own. The general
    // RDKit❗✔️:   # branch is kept because this function accepts an arbitrary matrix.
    // RDKit❗✔️:   if A.shape[0] == A.shape[1] and numpy.array_equal(A, A.T):
    // RDKit❗✔️:     roots = numpy.linalg.eigvalsh(A)
    // RDKit❗✔️:     res = numpy.array([1.0])
    // RDKit❗✔️:   else:
    // RDKit❗✔️:     roots = numpy.linalg.eigvals(A)
    // RDKit❗✔️:     res = numpy.array([1.0 + 0j])
    // RDKit❗✔️:   # Ordering the roots by increasing magnitude keeps every partial product near the scale of
    // RDKit❗✔️:   # the final coefficients, so nothing cancels. Multiplying them in arbitrary order (as
    // RDKit❗✔️:   # numpy.poly does) still reaches 4.8e4 relative error at 160 atoms.
    // RDKit❗✔️:   for r in roots[numpy.argsort(numpy.abs(roots))]:
    // RDKit❗✔️:     res = numpy.convolve(res, numpy.array([1.0, -r]))
    // RDKit❗✔️:   return numpy.real(res)
    // NumPy copies its LAPACK input to column-major storage once. Exact
    // symmetry is required: a one-triangle driver must never change an
    // asymmetric override. Both drivers own their detached work buffer.
    let n = graph.dimension;
    let mut symmetric = true;
    'rows: for i in 0..n {
        for j in 0..i {
            if graph.adjacency[i * n + j] != graph.adjacency[j * n + i] {
                symmetric = false;
                break 'rows;
            }
        }
    }
    let mut matrix = vec![0.0; graph.adjacency.len()];
    for j in 0..n {
        for i in 0..n {
            matrix[i + j * n] = graph.adjacency[i * n + j];
        }
    }
    if symmetric {
        let roots = lapack_symmetric::eigenvalues(n, matrix)?;
        let magnitudes: Vec<_> = roots.iter().map(|r| r.abs()).collect();
        let mut polynomial = vec![1.0];
        for index in numpy_argsort::argsort(&magnitudes) {
            polynomial = numpy_convolve::real_factor(&polynomial, roots[index]);
        }
        Ok(polynomial)
    } else {
        let roots = lapack_general::eigenvalues(n, matrix)?;
        // NumPy _linalg.py::_to_real_if_imag_zero:
        // NumPy❗✔️: if not isComplexType(result_t) and all(w.imag == 0.0):
        // NumPy❗✔️:     w = w.real
        // NumPy❗✔️:     result_t = _realType(result_t)
        // NumPy❗✔️: else:
        // NumPy❗✔️:     result_t = _complexType(result_t)
        // NumPy❗✔️: return w.astype(result_t, copy=False)
        let real_roots = roots.iter().all(|root| root[1] == 0.0);
        let magnitudes: Vec<_> = roots
            .iter()
            .copied()
            .map(|root| {
                if real_roots {
                    root[0].abs()
                } else {
                    numpy_convolve::complex_absolute(root)
                }
            })
            .collect();
        let mut polynomial = vec![[1.0, 0.0]];
        for index in numpy_argsort::argsort(&magnitudes) {
            polynomial = numpy_convolve::complex_factor(&polynomial, roots[index], real_roots);
        }
        Ok(polynomial.into_iter().map(|value| value[0]).collect())
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
    fn coefficients_follow_the_released_eigenvalue_expansion() {
        // RDKit 2026.09.2 Graphs.py / NumPy 2.5.3, fixed single-thread
        // OpenBLAS 446c436e: eigenvalue expansion retains rounding residues.
        let triangle = characteristic_polynomial(&graph(3, &[(0, 1), (1, 2), (0, 2)])).unwrap();
        assert_eq!(
            triangle.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
            [
                4607182418800017408,
                4377498837804122112,
                13837309855095848958,
                13835058055282163710
            ]
        );
        assert_eq!(
            characteristic_polynomial(&graph(4, &[(0, 1), (2, 3)])).unwrap(),
            [1.0, -0.0, -2.0, -0.0, 1.0]
        );
        assert_eq!(information_content(&graph(0, &[]), false).unwrap(), 0.0);
        assert_eq!(information_content(&graph(4, &[]), true).unwrap(), 0.0);
        let edges: Vec<_> = (0..119).map(|i| (i, i + 1)).collect();
        // The target eigenvalue expansion fixes the former long-chain
        // Faddeev instability without substituting ExactPolynomial.
        assert_eq!(
            information_content(&graph(120, &edges), true)
                .unwrap()
                .to_bits(),
            4615647539083575092
        );
    }
    #[test]
    fn asymmetric_overrides_use_the_whole_matrix() {
        let directed = graph_from_rows(
            2,
            [vec![0.0, 0.0], vec![1.0, 0.0]].iter().map(Vec::as_slice),
        )
        .unwrap();
        assert_eq!(information_content(&directed, false).unwrap(), 0.0);
        let looped = graph_from_rows(1, [vec![1.0]].iter().map(Vec::as_slice)).unwrap();
        assert_eq!(information_content(&looped, false).unwrap(), 2.0);
        assert_eq!(information_content(&looped, true).unwrap(), 1.0);
        // Directed 3-cycle has a complex-conjugate root pair. Keep the
        // complex polynomial until Graphs.py's final real projection.
        let cycle = graph_from_rows(
            3,
            [vec![0., 1., 0.], vec![0., 0., 1.], vec![1., 0., 0.]]
                .iter()
                .map(Vec::as_slice),
        )
        .unwrap();
        assert_eq!(
            characteristic_polynomial(&cycle)
                .unwrap()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            [
                4607182418800017408,
                4372995238176751616,
                4379750637617807360,
                13830554455654793217
            ]
        );
    }
}
