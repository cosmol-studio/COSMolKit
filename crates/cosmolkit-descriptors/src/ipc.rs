//! Exact characteristic-polynomial information content for detached molecular graphs.
//!
//! Modified port of the pinned RDKit implementation by COSMolKit contributors (2026).
//! cuik-molmaker portions: Copyright (c) 2025 NVIDIA CORPORATION & AFFILIATES.
//!
//! Numerical reconstruction follows cuik-molmaker 1c409b8fbcf428ca45a858ab62cf11035093c292
//! (`characteristic_polynomial.cpp`, Apache-2.0), with RDKit
//! b9518c8a8eb7ffc042e38e31d20d93bce33db10 (`CharacteristicPolynomial.cpp`, BSD-3-Clause)
//! defining loop handling, coefficient scaling, and Ipc/AvgIpc projection.
//! The applicable license notices are retained in the repository LICENSE.

use std::fmt;

use cosmolkit_model::TopologyBlock;
use num_bigint::BigUint;

mod entropy_log;

const LARGE_PRIMES: [u64; 64] = [
    0x3fffffffffffffc7,
    0x3fffffffffffffa9,
    0x3fffffffffffff8b,
    0x3fffffffffffff71,
    0x3fffffffffffff67,
    0x3fffffffffffff59,
    0x3fffffffffffff55,
    0x3fffffffffffff3d,
    0x3fffffffffffff35,
    0x3ffffffffffffeef,
    0x3ffffffffffffee1,
    0x3ffffffffffffec3,
    0x3ffffffffffffe45,
    0x3ffffffffffffe1d,
    0x3ffffffffffffe11,
    0x3ffffffffffffdc1,
    0x3ffffffffffffdbb,
    0x3ffffffffffffda5,
    0x3ffffffffffffd87,
    0x3ffffffffffffd69,
    0x3ffffffffffffd03,
    0x3ffffffffffffcfb,
    0x3ffffffffffffcf7,
    0x3ffffffffffffce9,
    0x3ffffffffffffcd3,
    0x3ffffffffffffcc1,
    0x3ffffffffffffc65,
    0x3ffffffffffffc2b,
    0x3ffffffffffffc1f,
    0x3ffffffffffffc17,
    0x3ffffffffffffc11,
    0x3ffffffffffffc07,
    0x3ffffffffffffb53,
    0x3ffffffffffffb27,
    0x3ffffffffffffaf3,
    0x3ffffffffffffab7,
    0x3ffffffffffffa67,
    0x3ffffffffffffa15,
    0x3ffffffffffff9ef,
    0x3ffffffffffff9d9,
    0x3ffffffffffff9d3,
    0x3ffffffffffff9c5,
    0x3ffffffffffff9af,
    0x3ffffffffffff977,
    0x3ffffffffffff95f,
    0x3ffffffffffff95b,
    0x3ffffffffffff959,
    0x3ffffffffffff8e1,
    0x3ffffffffffff8a7,
    0x3ffffffffffff889,
    0x3ffffffffffff87d,
    0x3ffffffffffff805,
    0x3ffffffffffff7e7,
    0x3ffffffffffff7c9,
    0x3ffffffffffff7a3,
    0x3ffffffffffff775,
    0x3ffffffffffff757,
    0x3ffffffffffff739,
    0x3ffffffffffff713,
    0x3ffffffffffff6d1,
    0x3ffffffffffff6c1,
    0x3ffffffffffff6b9,
    0x3ffffffffffff6a3,
    0x3ffffffffffff68b,
];

/// A supplied matrix is borrowed and interpreted from its lower triangle only.
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
    Ok(information_content(graph, false))
}

/// Computes total information content using a checked, borrowed matrix override.
pub fn ipc_with_distance_matrix(
    topology: &TopologyBlock,
    matrix: IpcMatrixInput<'_>,
) -> Result<f64, IpcError> {
    let graph = graph_from_matrix(topology.atoms.len(), matrix)?;
    Ok(information_content(graph, false))
}

/// Computes average information content from the supplied graph.
pub fn avg_ipc(topology: &TopologyBlock) -> Result<f64, IpcError> {
    let graph = graph_from_topology(topology)?;
    Ok(information_content(graph, true))
}

/// Computes average information content using a checked, borrowed matrix override.
pub fn avg_ipc_with_distance_matrix(
    topology: &TopologyBlock,
    matrix: IpcMatrixInput<'_>,
) -> Result<f64, IpcError> {
    let graph = graph_from_matrix(topology.atoms.len(), matrix)?;
    Ok(information_content(graph, true))
}

#[derive(Debug)]
struct BinaryGraph {
    num_vertices: usize,
    edges: Vec<(usize, usize)>,
    loops: Vec<bool>,
}

fn graph_from_topology(topology: &TopologyBlock) -> Result<BinaryGraph, IpcError> {
    let num_vertices = topology.atoms.len();
    let mut edges = Vec::with_capacity(topology.bonds.len());
    let mut loops = vec![false; num_vertices];

    for (bond_index, bond) in topology.bonds.iter().enumerate() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        for (endpoint, atom) in [("begin", begin), ("end", end)] {
            if atom >= num_vertices {
                return Err(IpcError::InvalidBondEndpoint {
                    bond: bond_index,
                    endpoint,
                    atom,
                    atom_count: num_vertices,
                });
            }
        }
        if begin == end {
            loops[begin] = true;
        } else {
            edges.push((begin, end));
        }
    }

    Ok(BinaryGraph {
        num_vertices,
        edges,
        loops,
    })
}

fn graph_from_matrix(
    num_vertices: usize,
    matrix: IpcMatrixInput<'_>,
) -> Result<BinaryGraph, IpcError> {
    match matrix {
        IpcMatrixInput::Rows(rows) => graph_from_rows(num_vertices, rows.iter().map(Vec::as_slice)),
        IpcMatrixInput::Dense(dense) => graph_from_rows(
            num_vertices,
            dense.values().chunks(dense.dimension().max(1)),
        ),
    }
}

fn graph_from_rows<'a>(
    num_vertices: usize,
    rows: impl ExactSizeIterator<Item = &'a [f64]>,
) -> Result<BinaryGraph, IpcError> {
    if rows.len() != num_vertices {
        return Err(IpcError::InvalidDistanceMatrixShape {
            row: None,
            actual: rows.len(),
            expected: num_vertices,
        });
    }

    let mut edges = Vec::new();
    let mut loops = Vec::with_capacity(num_vertices);
    for (row, values) in rows.enumerate() {
        if values.len() != num_vertices {
            return Err(IpcError::InvalidDistanceMatrixShape {
                row: Some(row),
                actual: values.len(),
                expected: num_vertices,
            });
        }
        loops.push(values[row] == 1.0);
        for (column, &value) in values[..row].iter().enumerate() {
            if value == 1.0 {
                edges.push((column, row));
            }
        }
    }

    Ok(BinaryGraph {
        num_vertices,
        edges,
        loops,
    })
}

fn information_content(graph: BinaryGraph, average: bool) -> f64 {
    // RDKit❗✔️: double calcIpc(const ROMol &mol, bool avg, const double *distMat) {
    // RDKit❗✔️:   const unsigned int nAts = mol.getNumAtoms();
    // RDKit❗✔️:   // The adjacency matrix: the bonds, or the entries of distMat equal to 1. Of
    // RDKit❗✔️:   // distMat only the lower triangle, diagonal included, is read, as the
    // RDKit❗✔️:   // symmetric eigensolver used before did.
    // RDKit❗✔️:   std::vector<std::pair<unsigned int, unsigned int>> edges;
    // RDKit❗✔️:   std::vector<bool> loops(nAts, false);
    // RDKit❗✔️:   if (distMat) {
    // RDKit❗✔️:     for (unsigned int i = 0; i < nAts; ++i) {
    // RDKit❗✔️:       for (unsigned int j = 0; j < i; ++j) {
    // RDKit❗✔️:         if (distMat[i * nAts + j] == 1.0) {
    // RDKit❗✔️:           edges.emplace_back(i, j);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       loops[i] = distMat[i * nAts + i] == 1.0;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:       edges.emplace_back(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // The coefficients are integers, computed exactly and rounded to double
    // RDKit❗✔️:   // once; see CharacteristicPolynomial.h.
    // RDKit❗✔️:   const auto cPoly =
    // RDKit❗✔️:       detail::characteristicPolynomialMagnitudes(nAts, edges, loops);
    // RDKit❗✔️:   if (!cPoly) {
    // RDKit❗✔️:     return std::numeric_limits<double>::quiet_NaN();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // The common power-of-two scale leaves the entropy unchanged; Ipc is
    // RDKit❗✔️:   // infinite once it leaves the double range.
    // RDKit❗✔️:   if (avg) {
    // RDKit❗✔️:     return RDInfoTheory::InfoEntropy(cPoly->scaled);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   double sum = 0.0;
    // RDKit❗✔️:   for (auto v : cPoly->scaled) {
    // RDKit❗✔️:     sum += v;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return std::ldexp(sum * RDInfoTheory::InfoEntropy(cPoly->scaled),
    // RDKit❗✔️:                     cPoly->exponent);
    // RDKit❗✔️: }
    let Some(magnitudes) = characteristic_polynomial_magnitudes(&graph) else {
        return f64::NAN;
    };
    let entropy = info_entropy(&magnitudes.scaled);
    if average {
        return entropy;
    }

    let mut sum = 0.0;
    for &value in &magnitudes.scaled {
        sum += value;
    }
    ldexp(sum * entropy, magnitudes.exponent)
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

fn ldexp(mut value: f64, mut exponent: i32) -> f64 {
    if value == 0.0 || !value.is_finite() {
        return value;
    }

    while exponent > 1023 {
        value *= f64::from_bits(2046_u64 << 52);
        exponent -= 1023;
    }
    while exponent < -1022 {
        // Leave 53 bits before the normal/subnormal boundary so that only the
        // final multiplication can round to a nonzero subnormal result.
        value *= f64::from_bits(54_u64 << 52); // 2^-969
        exponent += 969;
    }
    value * f64::from_bits(((exponent + 1023) as u64) << 52)
}

#[derive(Debug)]
struct PolynomialMagnitudes {
    scaled: Vec<f64>,
    exponent: i32,
}

#[derive(Debug)]
struct EliminationPlan {
    num_vertices: usize,
    num_entries: usize,
    edges: Vec<usize>,
    steps: Vec<EliminationStep>,
}

#[derive(Debug)]
struct EliminationStep {
    pivot: usize,
    to_later: Vec<usize>,
    later_diag: Vec<usize>,
    pairs: Vec<usize>,
}

#[derive(Debug, Clone, Copy)]
struct MontgomeryField {
    prime: u64,
    neg_inverse: u64,
    r_squared: u64,
    one: u64,
}

impl MontgomeryField {
    fn new(prime: u64) -> Self {
        // RDKit❗✔️: explicit MontgomeryField(std::uint64_t prime) : d_p(prime) {
        // RDKit❗✔️:   // Newton's iteration for p^-1 mod 2^64; each step doubles the valid bits
        // RDKit❗✔️:   std::uint64_t inverse = d_p;
        // RDKit❗✔️:   for (int i = 0; i < 6; ++i) {
        // RDKit❗✔️:     inverse *= 2 - d_p * inverse;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   d_negInverse = ~inverse + 1;
        // RDKit❗✔️:   // 2^64 mod p is the Montgomery form of 1; doubling it 64 more times gives
        // RDKit❗✔️:   // R^2 mod p, which converts into Montgomery form.
        // RDKit❗✔️:   d_one = (~d_p + 1) % d_p;
        // RDKit❗✔️:   d_rSquared = d_one;
        // RDKit❗✔️:   for (int i = 0; i < 64; ++i) {
        // RDKit❗✔️:     d_rSquared = add(d_rSquared, d_rSquared);
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        let mut inverse = prime;
        for _ in 0..6 {
            inverse = inverse.wrapping_mul(2_u64.wrapping_sub(prime.wrapping_mul(inverse)));
        }
        let neg_inverse = inverse.wrapping_neg();
        let r = ((1_u128 << 64) % u128::from(prime)) as u64;
        let r_squared = ((u128::from(r) * u128::from(r)) % u128::from(prime)) as u64;
        let mut field = Self {
            prime,
            neg_inverse,
            r_squared,
            one: 0,
        };
        field.one = field.to(1);
        field
    }

    fn prime(self) -> u64 {
        self.prime
    }

    fn one(self) -> u64 {
        self.one
    }

    fn to(self, value: u64) -> u64 {
        // RDKit❗✔️: std::uint64_t to(std::uint64_t a) const { return mul(a % d_p, d_rSquared); }
        self.mul(value % self.prime, self.r_squared)
    }

    fn from(self, value: u64) -> u64 {
        // RDKit❗✔️: std::uint64_t from(std::uint64_t a) const { return reduce(0, a); }
        self.reduce(u128::from(value))
    }

    fn mul(self, left: u64, right: u64) -> u64 {
        // RDKit❗✔️: std::uint64_t mul(std::uint64_t a, std::uint64_t b) const {
        // RDKit❗✔️:   std::uint64_t high, low;
        // RDKit❗✔️:   multiplyWide(a, b, high, low);
        // RDKit❗✔️:   return reduce(high, low);
        // RDKit❗✔️: }
        self.reduce(u128::from(left) * u128::from(right))
    }

    fn add(self, left: u64, right: u64) -> u64 {
        // RDKit❗✔️: std::uint64_t add(std::uint64_t a, std::uint64_t b) const {
        // RDKit❗✔️:   const std::uint64_t s = a + b;  // < 2^63, no overflow
        // RDKit❗✔️:   return s >= d_p ? s - d_p : s;
        // RDKit❗✔️: }
        let sum = left + right;
        if sum >= self.prime {
            sum - self.prime
        } else {
            sum
        }
    }

    fn sub(self, left: u64, right: u64) -> u64 {
        // RDKit❗✔️: std::uint64_t sub(std::uint64_t a, std::uint64_t b) const {
        // RDKit❗✔️:   return a >= b ? a - b : a + d_p - b;
        // RDKit❗✔️: }
        if left >= right {
            left - right
        } else {
            left + self.prime - right
        }
    }

    fn inverse(self, mut value: u64) -> u64 {
        // RDKit❗✔️: std::uint64_t inverse(std::uint64_t a) const {
        // RDKit❗✔️:   std::uint64_t result = d_one;
        // RDKit❗✔️:   for (std::uint64_t e = d_p - 2; e; e >>= 1, a = mul(a, a)) {
        // RDKit❗✔️:     if (e & 1) {
        // RDKit❗✔️:       result = mul(result, a);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return result;
        // RDKit❗✔️: }
        let mut result = self.one;
        let mut exponent = self.prime - 2;
        while exponent != 0 {
            if exponent & 1 != 0 {
                result = self.mul(result, value);
            }
            value = self.mul(value, value);
            exponent >>= 1;
        }
        result
    }

    fn reduce(self, value: u128) -> u64 {
        // RDKit❗✔️: std::uint64_t reduce(std::uint64_t high, std::uint64_t low) const {
        // RDKit❗✔️:   const std::uint64_t m = low * d_negInverse;
        // RDKit❗✔️:   std::uint64_t mpHigh, mpLow;
        // RDKit❗✔️:   multiplyWide(m, d_p, mpHigh, mpLow);
        // RDKit❗✔️:   const std::uint64_t r = high + mpHigh + (low != 0 ? 1 : 0);
        // RDKit❗✔️:   return r >= d_p ? r - d_p : r;
        // RDKit❗✔️: }
        let m = (value as u64).wrapping_mul(self.neg_inverse);
        let reduced = ((value + u128::from(m) * u128::from(self.prime)) >> 64) as u64;
        if reduced >= self.prime {
            reduced - self.prime
        } else {
            reduced
        }
    }
}

fn plan_elimination(num_vertices: usize, edges: &[(usize, usize)]) -> EliminationPlan {
    // RDKit❗✔️: EliminationPlan planElimination(
    // RDKit❗✔️:     unsigned int numVertices,
    // RDKit❗✔️:     const std::vector<std::pair<unsigned int, unsigned int>> &edges) {
    // RDKit❗✔️:   EliminationPlan plan;
    // RDKit❗✔️:   plan.numVertices = numVertices;
    // RDKit❗✔️:   plan.numEntries = numVertices;
    // RDKit❗✔️:   std::vector<std::vector<std::pair<size_t, size_t>>> adjacency(numVertices);
    // RDKit❗✔️:   const auto connect = [&](size_t a, size_t b) {
    // RDKit❗✔️:     for (const auto &[vertex, entry] : adjacency[a]) {
    // RDKit❗✔️:       if (vertex == b) {
    // RDKit❗✔️:         return entry;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     const size_t entry = plan.numEntries++;
    // RDKit❗✔️:     adjacency[a].emplace_back(b, entry);
    // RDKit❗✔️:     adjacency[b].emplace_back(a, entry);
    // RDKit❗✔️:     return entry;
    // RDKit❗✔️:   };
    // RDKit❗✔️:   for (const auto &[a, b] : edges) {
    // RDKit❗✔️:     plan.edges.push_back(connect(a, b));
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<bool> eliminated(numVertices, false);
    // RDKit❗✔️:   for (size_t step = 0; step < numVertices; ++step) {
    // RDKit❗✔️:     size_t v = numVertices;
    // RDKit❗✔️:     for (size_t u = 0; u < numVertices; ++u) {
    // RDKit❗✔️:       if (!eliminated[u] &&
    // RDKit❗✔️:           (v == numVertices || adjacency[u].size() < adjacency[v].size())) {
    // RDKit❗✔️:         v = u;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     EliminationPlan::Step current{v, {}, {}, {}};
    // RDKit❗✔️:     const auto later = adjacency[v];
    // RDKit❗✔️:     for (const auto &[vertex, entry] : later) {
    // RDKit❗✔️:       current.toLater.push_back(entry);
    // RDKit❗✔️:       current.laterDiag.push_back(vertex);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (size_t i = 0; i < later.size(); ++i) {
    // RDKit❗✔️:       for (size_t j = i + 1; j < later.size(); ++j) {
    // RDKit❗✔️:         current.pairs.push_back(connect(later[i].first, later[j].first));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (const auto &[vertex, entry] : later) {
    // RDKit❗✔️:       auto &list = adjacency[vertex];
    // RDKit❗✔️:       list.erase(std::find_if(list.begin(), list.end(), [v](const auto &item) {
    // RDKit❗✔️:         return item.first == v;
    // RDKit❗✔️:       }));
    // RDKit❗✔️:     }
    // RDKit❗✔️:     adjacency[v].clear();
    // RDKit❗✔️:     eliminated[v] = true;
    // RDKit❗✔️:     plan.steps.push_back(std::move(current));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return plan;
    // RDKit❗✔️: }
    let mut plan = EliminationPlan {
        num_vertices,
        num_entries: num_vertices,
        edges: Vec::with_capacity(edges.len()),
        steps: Vec::with_capacity(num_vertices),
    };
    let mut adjacency = vec![Vec::<(usize, usize)>::new(); num_vertices];

    fn connect(
        adjacency: &mut [Vec<(usize, usize)>],
        entries: &mut usize,
        a: usize,
        b: usize,
    ) -> usize {
        if let Some((_, entry)) = adjacency[a].iter().find(|(vertex, _)| *vertex == b) {
            return *entry;
        }
        let entry = *entries;
        *entries += 1;
        adjacency[a].push((b, entry));
        adjacency[b].push((a, entry));
        entry
    }

    for &(a, b) in edges {
        let entry = connect(&mut adjacency, &mut plan.num_entries, a, b);
        plan.edges.push(entry);
    }

    let mut eliminated = vec![false; num_vertices];
    for _ in 0..num_vertices {
        let mut vertex = num_vertices;
        for candidate in 0..num_vertices {
            if !eliminated[candidate]
                && (vertex == num_vertices || adjacency[candidate].len() < adjacency[vertex].len())
            {
                vertex = candidate;
            }
        }

        let later = adjacency[vertex].clone();
        let mut step = EliminationStep {
            pivot: vertex,
            to_later: Vec::with_capacity(later.len()),
            later_diag: Vec::with_capacity(later.len()),
            pairs: Vec::new(),
        };
        for &(later_vertex, entry) in &later {
            step.to_later.push(entry);
            step.later_diag.push(later_vertex);
        }
        for left in 0..later.len() {
            for right in left + 1..later.len() {
                let entry = connect(
                    &mut adjacency,
                    &mut plan.num_entries,
                    later[left].0,
                    later[right].0,
                );
                step.pairs.push(entry);
            }
        }
        for &(later_vertex, _) in &later {
            adjacency[later_vertex].retain(|(neighbor, _)| *neighbor != vertex);
        }
        adjacency[vertex].clear();
        eliminated[vertex] = true;
        plan.steps.push(step);
    }
    plan
}

fn evaluate_determinants(
    plan: &EliminationPlan,
    loops: &[bool],
    field: MontgomeryField,
    first: u64,
    count: usize,
) -> (Vec<u64>, Vec<bool>) {
    // RDKit❗✔️: void evaluateDeterminants(const EliminationPlan &plan,
    // RDKit❗✔️:                           const std::vector<bool> &loops,
    // RDKit❗✔️:                           const MontgomeryField &field, std::uint64_t first,
    // RDKit❗✔️:                           size_t count,
    // RDKit❗✔️:                           std::vector<std::uint64_t> &determinants,
    // RDKit❗✔️:                           std::vector<bool> &usable) {
    // RDKit❗✔️:   std::vector<std::uint64_t> values(plan.numEntries * count, 0);
    // RDKit❗✔️:   std::vector<std::uint64_t> inverse(count), prefix(count), multipliers;
    // RDKit❗✔️:   for (size_t v = 0; v < plan.numVertices; ++v) {
    // RDKit❗✔️:     for (size_t j = 0; j < count; ++j) {
    // RDKit❗✔️:       const std::uint64_t x = field.to(first + j);
    // RDKit❗✔️:       values[v * count + j] = loops[v] ? field.sub(x, field.one()) : x;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const std::uint64_t minusOne = field.sub(0, field.one());
    // RDKit❗✔️:   for (const size_t edge : plan.edges) {
    // RDKit❗✔️:     std::fill_n(&values[edge * count], count, minusOne);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   determinants.assign(count, field.one());
    // RDKit❗✔️:   usable.assign(count, true);
    // RDKit❗✔️:
    // RDKit❗✔️:   for (const auto &step : plan.steps) {
    // RDKit❗✔️:     std::uint64_t *const pivot = &values[step.pivot * count];
    // RDKit❗✔️:     for (size_t j = 0; j < count; ++j) {
    // RDKit❗✔️:       if (pivot[j] == 0) {
    // RDKit❗✔️:         usable[j] = false;
    // RDKit❗✔️:         pivot[j] = field.one();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     prefix[0] = pivot[0];
    // RDKit❗✔️:     for (size_t j = 1; j < count; ++j) {
    // RDKit❗✔️:       prefix[j] = field.mul(prefix[j - 1], pivot[j]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     std::uint64_t running = field.inverse(prefix[count - 1]);
    // RDKit❗✔️:     for (size_t j = count; j-- > 1;) {
    // RDKit❗✔️:       inverse[j] = field.mul(running, prefix[j - 1]);
    // RDKit❗✔️:       running = field.mul(running, pivot[j]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     inverse[0] = running;
    // RDKit❗✔️:     for (size_t j = 0; j < count; ++j) {
    // RDKit❗✔️:       determinants[j] = field.mul(determinants[j], pivot[j]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     const size_t later = step.toLater.size();
    // RDKit❗✔️:     multipliers.assign(later * count, 0);
    // RDKit❗✔️:     for (size_t a = 0; a < later; ++a) {
    // RDKit❗✔️:       const std::uint64_t *const column = &values[step.toLater[a] * count];
    // RDKit❗✔️:       std::uint64_t *const factor = &multipliers[a * count];
    // RDKit❗✔️:       std::uint64_t *const diagonal = &values[step.laterDiag[a] * count];
    // RDKit❗✔️:       for (size_t j = 0; j < count; ++j) {
    // RDKit❗✔️:         factor[j] = field.mul(column[j], inverse[j]);
    // RDKit❗✔️:         diagonal[j] = field.sub(diagonal[j], field.mul(factor[j], column[j]));
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     size_t pair = 0;
    // RDKit❗✔️:     for (size_t a = 0; a < later; ++a) {
    // RDKit❗✔️:       for (size_t b = a + 1; b < later; ++b, ++pair) {
    // RDKit❗✔️:         const std::uint64_t *const factor = &multipliers[a * count];
    // RDKit❗✔️:         const std::uint64_t *const column = &values[step.toLater[b] * count];
    // RDKit❗✔️:         std::uint64_t *const target = &values[step.pairs[pair] * count];
    // RDKit❗✔️:         for (size_t j = 0; j < count; ++j) {
    // RDKit❗✔️:           target[j] = field.sub(target[j], field.mul(factor[j], column[j]));
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    let mut values = vec![0_u64; plan.num_entries * count];
    let mut inverse = vec![0_u64; count];
    let mut prefix = vec![0_u64; count];
    let mut multipliers = Vec::new();

    for vertex in 0..plan.num_vertices {
        for point in 0..count {
            let x = field.to(first + point as u64);
            values[vertex * count + point] = if loops[vertex] {
                field.sub(x, field.one())
            } else {
                x
            };
        }
    }
    let minus_one = field.sub(0, field.one());
    for &edge in &plan.edges {
        values[edge * count..(edge + 1) * count].fill(minus_one);
    }

    let mut determinants = vec![field.one(); count];
    let mut usable = vec![true; count];
    for step in &plan.steps {
        let pivot_start = step.pivot * count;
        for point in 0..count {
            if values[pivot_start + point] == 0 {
                usable[point] = false;
                values[pivot_start + point] = field.one();
            }
        }
        prefix[0] = values[pivot_start];
        for point in 1..count {
            prefix[point] = field.mul(prefix[point - 1], values[pivot_start + point]);
        }
        let mut running = field.inverse(prefix[count - 1]);
        for point in (1..count).rev() {
            inverse[point] = field.mul(running, prefix[point - 1]);
            running = field.mul(running, values[pivot_start + point]);
        }
        inverse[0] = running;
        for point in 0..count {
            determinants[point] = field.mul(determinants[point], values[pivot_start + point]);
        }

        let later = step.to_later.len();
        multipliers.resize(later * count, 0);
        for neighbor in 0..later {
            let column_start = step.to_later[neighbor] * count;
            let factor_start = neighbor * count;
            let diagonal_start = step.later_diag[neighbor] * count;
            for point in 0..count {
                let factor = field.mul(values[column_start + point], inverse[point]);
                multipliers[factor_start + point] = factor;
                values[diagonal_start + point] = field.sub(
                    values[diagonal_start + point],
                    field.mul(factor, values[column_start + point]),
                );
            }
        }
        let mut pair = 0;
        for left in 0..later {
            for right in left + 1..later {
                let factor_start = left * count;
                let column_start = step.to_later[right] * count;
                let target_start = step.pairs[pair] * count;
                for point in 0..count {
                    values[target_start + point] = field.sub(
                        values[target_start + point],
                        field.mul(
                            multipliers[factor_start + point],
                            values[column_start + point],
                        ),
                    );
                }
                pair += 1;
            }
        }
    }

    (determinants, usable)
}

fn coefficients_mod_prime(
    plan: &EliminationPlan,
    loops: &[bool],
    field: MontgomeryField,
    first_point: u64,
) -> Vec<u64> {
    // RDKit❗✔️: std::vector<std::uint64_t> coefficientsModPrime(const EliminationPlan &plan,
    // RDKit❗✔️:                                                 const std::vector<bool> &loops,
    // RDKit❗✔️:                                                 const MontgomeryField &field,
    // RDKit❗✔️:                                                 std::uint64_t firstPoint) {
    // RDKit❗✔️:   const size_t n = plan.numVertices;
    // RDKit❗✔️:   constexpr size_t maxBlock = 256;
    // RDKit❗✔️:   std::vector<std::uint64_t> xs, ys, determinants;
    // RDKit❗✔️:   std::vector<bool> usable;
    // RDKit❗✔️:   for (std::uint64_t next = firstPoint; xs.size() < n + 1;) {
    // RDKit❗✔️:     const size_t count = std::min<size_t>(maxBlock, n + 1 - xs.size() + 4);
    // RDKit❗✔️:     evaluateDeterminants(plan, loops, field, next, count, determinants, usable);
    // RDKit❗✔️:     for (size_t j = 0; j < count && xs.size() < n + 1; ++j) {
    // RDKit❗✔️:       if (usable[j]) {
    // RDKit❗✔️:         xs.push_back(next + j);
    // RDKit❗✔️:         ys.push_back(determinants[j]);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     next += count;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   const std::uint64_t p = field.prime();
    // RDKit❗✔️:   const size_t span = static_cast<size_t>(xs.back() - xs.front());
    // RDKit❗✔️:   std::vector<std::uint64_t> smallInverse(span + 1, 0);
    // RDKit❗✔️:   if (span >= 1) {
    // RDKit❗✔️:     smallInverse[1] = field.one();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (size_t i = 2; i <= span; ++i) {
    // RDKit❗✔️:     smallInverse[i] =
    // RDKit❗✔️:         field.sub(0, field.mul(field.to(p / i), smallInverse[p % i]));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (size_t k = 1; k <= n; ++k) {
    // RDKit❗✔️:     for (size_t i = n; i >= k; --i) {
    // RDKit❗✔️:       ys[i] = field.mul(field.sub(ys[i], ys[i - 1]),
    // RDKit❗✔️:                         smallInverse[xs[i] - xs[i - k]]);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<std::uint64_t> polynomial(n + 1, 0);
    // RDKit❗✔️:   polynomial[0] = ys[n];
    // RDKit❗✔️:   for (size_t i = n, degree = 0; i-- > 0; ++degree) {
    // RDKit❗✔️:     const std::uint64_t xi = field.to(xs[i]);
    // RDKit❗✔️:     for (size_t t = degree + 1; t > 0; --t) {
    // RDKit❗✔️:       polynomial[t] =
    // RDKit❗✔️:           field.sub(polynomial[t - 1], field.mul(xi, polynomial[t]));
    // RDKit❗✔️:     }
    // RDKit❗✔️:     polynomial[0] =
    // RDKit❗✔️:         field.add(field.sub(0, field.mul(xi, polynomial[0])), ys[i]);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto &c : polynomial) {
    // RDKit❗✔️:     c = field.from(c);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return polynomial;
    // RDKit❗✔️: }
    let n = plan.num_vertices;
    let mut xs = Vec::with_capacity(n + 1);
    let mut ys = Vec::with_capacity(n + 1);
    let mut next = first_point;

    while xs.len() < n + 1 {
        let count = 256.min(n + 1 - xs.len() + 4);
        let (determinants, usable) = evaluate_determinants(plan, loops, field, next, count);
        for point in 0..count {
            if usable[point] && xs.len() < n + 1 {
                xs.push(next + point as u64);
                ys.push(determinants[point]);
            }
        }
        next += count as u64;
    }

    let span = (xs[n] - xs[0]) as usize;
    let prime = field.prime();
    let mut small_inverse = vec![0_u64; span + 1];
    if span >= 1 {
        small_inverse[1] = field.one();
    }
    for gap in 2..=span {
        let quotient = prime / gap as u64;
        let remainder = (prime % gap as u64) as usize;
        small_inverse[gap] = field.sub(0, field.mul(field.to(quotient), small_inverse[remainder]));
    }

    for order in 1..=n {
        for i in (order..=n).rev() {
            let denominator = (xs[i] - xs[i - order]) as usize;
            ys[i] = field.mul(field.sub(ys[i], ys[i - 1]), small_inverse[denominator]);
        }
    }

    let mut polynomial = vec![0_u64; n + 1];
    polynomial[0] = ys[n];
    for i in (0..n).rev() {
        let degree = n - i - 1;
        let x = field.to(xs[i]);
        for term in (1..=degree + 1).rev() {
            polynomial[term] = field.sub(polynomial[term - 1], field.mul(x, polynomial[term]));
        }
        polynomial[0] = field.add(field.sub(0, field.mul(x, polynomial[0])), ys[i]);
    }
    for coefficient in &mut polynomial {
        *coefficient = field.from(*coefficient);
    }
    polynomial
}

fn coefficient_bits_bound(num_vertices: usize, trace_of_square: usize) -> f64 {
    // RDKit❗✔️: double coefficientBitsBound(size_t n, size_t traceOfSquare) {
    // RDKit❗✔️:   const double halfLogDensity =
    // RDKit❗✔️:       0.5 * std::log2(static_cast<double>(traceOfSquare) / n);
    // RDKit❗✔️:   const double log2NFactorial = std::lgamma(n + 1.0) / std::log(2.0);
    // RDKit❗✔️:   double bound = 0.0;
    // RDKit❗✔️:   for (size_t k = 0; k <= n; ++k) {
    // RDKit❗✔️:     const double log2Binomial =
    // RDKit❗✔️:         log2NFactorial -
    // RDKit❗✔️:         (std::lgamma(k + 1.0) + std::lgamma(n - k + 1.0)) / std::log(2.0);
    // RDKit❗✔️:     bound = std::max(bound, log2Binomial + k * halfLogDensity);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return bound;
    // RDKit❗✔️: }
    let n = num_vertices as f64;
    let half_log_density = 0.5 * (trace_of_square as f64 / n).log2();
    let log2_n_factorial = log_gamma(num_vertices as f64 + 1.0) / std::f64::consts::LN_2;
    let mut bound: f64 = 0.0;
    for k in 0..=num_vertices {
        let log2_binomial = log2_n_factorial
            - (log_gamma(k as f64 + 1.0) + log_gamma((num_vertices - k) as f64 + 1.0))
                / std::f64::consts::LN_2;
        bound = bound.max(log2_binomial + k as f64 * half_log_density);
    }
    bound
}

fn log_gamma(value: f64) -> f64 {
    // Lanczos at positive integer arguments, used only to select the prime count.
    const COEFFICIENTS: [f64; 8] = [
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    let z = value - 1.0;
    let mut sum = 0.9999999999998099;
    for (index, coefficient) in COEFFICIENTS.iter().enumerate() {
        sum += coefficient / (z + index as f64 + 1.0);
    }
    let t = z + COEFFICIENTS.len() as f64 - 0.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + sum.ln()
}

fn magnitude(value: &BigUint, exponent: i32) -> f64 {
    // RDKit❗✔️: double magnitude(BigInt value, int exponent) {
    // RDKit❗✔️:   if (value < 0) {
    // RDKit❗✔️:     value = -value;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (value == 0) {
    // RDKit❗✔️:     return 0.0;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const unsigned int top = boost::multiprecision::msb(value);
    // RDKit❗✔️:   if (top < 64) {
    // RDKit❗✔️:     return std::ldexp(static_cast<double>(static_cast<std::uint64_t>(value)),
    // RDKit❗✔️:                       -exponent);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const unsigned int shift = top - 63;
    // RDKit❗✔️:   auto head = static_cast<std::uint64_t>(value >> shift);
    // RDKit❗✔️:   if (boost::multiprecision::lsb(value) < shift) {
    // RDKit❗✔️:     head |= 1u;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return std::ldexp(static_cast<double>(head),
    // RDKit❗✔️:                     static_cast<int>(shift) - exponent);
    // RDKit❗✔️: }
    let Some(low) = value.trailing_zeros() else {
        return 0.0;
    };
    let top = value.bits() as usize - 1;
    if top < 64 {
        let head = value.iter_u64_digits().next().unwrap_or(0);
        return ldexp(head as f64, -exponent);
    }

    let shift = top - 63;
    let mut head = (value >> shift).iter_u64_digits().next().unwrap_or(0);
    if low < shift as u64 {
        head |= 1;
    }
    ldexp(head as f64, shift as i32 - exponent)
}

fn characteristic_polynomial_magnitudes(graph: &BinaryGraph) -> Option<PolynomialMagnitudes> {
    // RDKit❗✔️: std::optional<PolynomialMagnitudes> characteristicPolynomialMagnitudes(
    // RDKit❗✔️:     unsigned int numVertices,
    // RDKit❗✔️:     const std::vector<std::pair<unsigned int, unsigned int>> &edges,
    // RDKit❗✔️:     const std::vector<bool> &loops) {
    // RDKit❗✔️:   const size_t n = numVertices;
    // RDKit❗✔️:   PolynomialMagnitudes result{std::vector<double>(n + 1, 0.0)};
    // RDKit❗✔️:   result.scaled[0] = 1.0;
    // RDKit❗✔️:   std::vector<size_t> rowSums(n, 0);
    // RDKit❗✔️:   size_t traceOfSquare = 0;  // tr A^2 = 2 |edges| + |loops|
    // RDKit❗✔️:   for (const auto &[a, b] : edges) {
    // RDKit❗✔️:     ++rowSums[a];
    // RDKit❗✔️:     ++rowSums[b];
    // RDKit❗✔️:     traceOfSquare += 2;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (size_t v = 0; v < n; ++v) {
    // RDKit❗✔️:     if (loops[v]) {
    // RDKit❗✔️:       ++rowSums[v];
    // RDKit❗✔️:       ++traceOfSquare;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (traceOfSquare == 0) {
    // RDKit❗✔️:     return result;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   const double neededBits = coefficientBitsBound(n, traceOfSquare) + 3.0;
    // RDKit❗✔️:   size_t numPrimes = 0;
    // RDKit❗✔️:   for (double covered = 0.0; covered < neededBits; ++numPrimes) {
    // RDKit❗✔️:     if (numPrimes == largePrimes.size()) {
    // RDKit❗✔️:       return std::nullopt;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     covered += std::log2(static_cast<double>(largePrimes[numPrimes]));
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   const EliminationPlan plan = planElimination(numVertices, edges);
    // RDKit❗✔️:   const std::uint64_t first =
    // RDKit❗✔️:       *std::max_element(rowSums.begin(), rowSums.end()) + 1;
    // RDKit❗✔️:   const auto residues =
    // RDKit❗✔️:       coefficientsModPrime(plan, loops, MontgomeryField(largePrimes[0]), first);
    // RDKit❗✔️:
    // RDKit❗✔️:   if (numPrimes == 1) {
    // RDKit❗✔️:     const std::uint64_t p = largePrimes[0];
    // RDKit❗✔️:     for (size_t k = 0; k <= n; ++k) {
    // RDKit❗✔️:       const std::uint64_t r = residues[n - k];
    // RDKit❗✔️:       result.scaled[k] = static_cast<double>(r > p / 2 ? p - r : r);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     return result;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<BigInt> exact(residues.begin(), residues.end());
    // RDKit❗✔️:   BigInt modulus = largePrimes[0];
    // RDKit❗✔️:   for (size_t i = 1; i < numPrimes; ++i) {
    // RDKit❗✔️:     const std::uint64_t p = largePrimes[i];
    // RDKit❗✔️:     const MontgomeryField field(p);
    // RDKit❗✔️:     const auto r = coefficientsModPrime(plan, loops, field, first);
    // RDKit❗✔️:     const std::uint64_t inverseOfModulus =
    // RDKit❗✔️:         field.inverse(field.to(static_cast<std::uint64_t>(modulus % p)));
    // RDKit❗✔️:     for (size_t k = 0; k <= n; ++k) {
    // RDKit❗✔️:       const auto current = static_cast<std::uint64_t>(exact[k] % p);
    // RDKit❗✔️:       const std::uint64_t delta =
    // RDKit❗✔️:           r[k] >= current ? r[k] - current : r[k] + p - current;
    // RDKit❗✔️:       exact[k] +=
    // RDKit❗✔️:           modulus * field.from(field.mul(field.to(delta), inverseOfModulus));
    // RDKit❗✔️:     }
    // RDKit❗✔️:     modulus *= p;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   const BigInt half = modulus / 2;
    // RDKit❗✔️:   unsigned int top = 0;
    // RDKit❗✔️:   for (BigInt &x : exact) {
    // RDKit❗✔️:     if (x > half) {
    // RDKit❗✔️:       x = modulus - x;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (x != 0) {
    // RDKit❗✔️:       top = std::max(top, static_cast<unsigned int>(boost::multiprecision::msb(x)));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const int sumBits =
    // RDKit❗✔️:       static_cast<int>(top) + 2 + static_cast<int>(std::bit_width(n + 1));
    // RDKit❗✔️:   result.exponent =
    // RDKit❗✔️:       std::max(0, sumBits - std::numeric_limits<double>::max_exponent);
    // RDKit❗✔️:   for (size_t k = 0; k <= n; ++k) {
    // RDKit❗✔️:     result.scaled[k] = magnitude(exact[n - k], result.exponent);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return result;
    // RDKit❗✔️: }
    let n = graph.num_vertices;
    let mut result = PolynomialMagnitudes {
        scaled: vec![0.0; n + 1],
        exponent: 0,
    };
    result.scaled[0] = 1.0;
    if n == 0 {
        return Some(result);
    }

    let mut row_sums = vec![0_usize; n];
    let mut trace_of_square = 0_usize;
    for &(a, b) in &graph.edges {
        row_sums[a] += 1;
        row_sums[b] += 1;
        trace_of_square += 2;
    }
    for (vertex, has_loop) in graph.loops.iter().copied().enumerate() {
        if has_loop {
            row_sums[vertex] += 1;
            trace_of_square += 1;
        }
    }
    if trace_of_square == 0 {
        return Some(result);
    }

    let needed_bits = coefficient_bits_bound(n, trace_of_square) + 3.0;
    let mut covered = 0.0;
    let mut num_primes = 0;
    while covered < needed_bits {
        let Some(&prime) = LARGE_PRIMES.get(num_primes) else {
            return None;
        };
        covered += (prime as f64).log2();
        num_primes += 1;
    }

    let plan = plan_elimination(n, &graph.edges);
    let first = row_sums.iter().copied().max().unwrap_or(0) as u64 + 1;
    let residues = coefficients_mod_prime(
        &plan,
        &graph.loops,
        MontgomeryField::new(LARGE_PRIMES[0]),
        first,
    );
    if num_primes == 1 {
        let prime = LARGE_PRIMES[0];
        for k in 0..=n {
            let residue = residues[n - k];
            result.scaled[k] = (if residue > prime / 2 {
                prime - residue
            } else {
                residue
            }) as f64;
        }
        return Some(result);
    }

    let mut exact = residues
        .iter()
        .copied()
        .map(BigUint::from)
        .collect::<Vec<_>>();
    let mut modulus = BigUint::from(LARGE_PRIMES[0]);
    for &prime in &LARGE_PRIMES[1..num_primes] {
        let field = MontgomeryField::new(prime);
        let modulus_mod_prime = (&modulus % prime).iter_u64_digits().next().unwrap_or(0);
        let inverse_of_modulus = field.inverse(field.to(modulus_mod_prime));
        let residues = coefficients_mod_prime(&plan, &graph.loops, field, first);
        for k in 0..=n {
            let current = (&exact[k] % prime).iter_u64_digits().next().unwrap_or(0);
            let delta = if residues[k] >= current {
                residues[k] - current
            } else {
                residues[k] + prime - current
            };
            let term = field.from(field.mul(field.to(delta), inverse_of_modulus));
            exact[k] += &modulus * term;
        }
        modulus *= prime;
    }

    let half = &modulus >> 1;
    let mut top = 0_usize;
    for coefficient in &mut exact {
        if *coefficient > half {
            *coefficient = &modulus - &*coefficient;
        }
        if coefficient != &BigUint::default() {
            top = top.max(coefficient.bits() as usize - 1);
        }
    }
    let sum_bits = top as i32 + 2 + (usize::BITS - n.saturating_add(1).leading_zeros()) as i32;
    result.exponent = (sum_bits - 1024).max(0);
    for k in 0..=n {
        result.scaled[k] = magnitude(&exact[n - k], result.exponent);
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graph(num_vertices: usize, edges: &[(usize, usize)]) -> BinaryGraph {
        BinaryGraph {
            num_vertices,
            edges: edges.to_vec(),
            loops: vec![false; num_vertices],
        }
    }

    fn path_graph(num_vertices: usize) -> BinaryGraph {
        graph(
            num_vertices,
            &(0..num_vertices.saturating_sub(1))
                .map(|i| (i, i + 1))
                .collect::<Vec<_>>(),
        )
    }

    #[test]
    fn coefficients_match_triangle_and_disconnected_edges() {
        let field = MontgomeryField::new(LARGE_PRIMES[0]);
        assert_eq!(field.from(field.one()), 1);
        assert_eq!(field.from(field.to(7)), 7);
        assert_eq!(field.mul(7, field.one()), 7);
        assert_eq!(field.mul(field.to(7), 1), 7);
        let triangle_graph = graph(3, &[(0, 1), (1, 2), (0, 2)]);
        let plan = plan_elimination(3, &triangle_graph.edges);
        let (determinants, usable) =
            evaluate_determinants(&plan, &triangle_graph.loops, field, 3, 8);
        assert_eq!(usable, [true; 8]);
        assert_eq!(
            determinants
                .iter()
                .map(|&value| field.from(value))
                .collect::<Vec<_>>(),
            [16, 50, 108, 196, 320, 486, 700, 968]
        );
        assert_eq!(
            coefficients_mod_prime(&plan, &triangle_graph.loops, field, 3),
            [LARGE_PRIMES[0] - 2, LARGE_PRIMES[0] - 3, 0, 1]
        );
        let triangle = characteristic_polynomial_magnitudes(&triangle_graph).unwrap();
        assert_eq!(triangle.scaled, [1.0, 0.0, 3.0, 2.0]);

        let disconnected =
            characteristic_polynomial_magnitudes(&graph(4, &[(0, 1), (2, 3)])).unwrap();
        assert_eq!(disconnected.scaled, [1.0, 0.0, 2.0, 0.0, 1.0]);
    }

    #[test]
    fn path_120_matches_independent_binomial_coefficients() {
        let actual = characteristic_polynomial_magnitudes(&path_graph(120)).unwrap();
        let mut expected = vec![0.0; 121];
        expected[0] = 1.0;
        for k in 1..=60 {
            let coefficient = binomial(120 - k, k) as f64;
            expected[2 * k] = coefficient;
        }
        assert_eq!(actual.scaled, expected);
    }

    #[test]
    fn empty_edgeless_and_unsupported_graphs_have_source_outcomes() {
        assert_eq!(information_content(graph(0, &[]), false), 0.0);
        assert_eq!(information_content(graph(0, &[]), true), 0.0);
        assert_eq!(information_content(graph(4_000, &[]), false), 0.0);
        assert_eq!(information_content(graph(4_000, &[]), true), 0.0);
        assert!(information_content(path_graph(4_000), false).is_nan());
        assert!(information_content(path_graph(4_000), true).is_nan());
    }

    #[test]
    fn binary_scaling_preserves_representable_results() {
        assert_eq!(ldexp(0.25, 1024), 2.0_f64.powi(1022));
        assert_eq!(ldexp(0.0, 2048), 0.0);
        assert_eq!(ldexp(2.0_f64.powi(1023), -2097).to_bits(), 1);
        assert_eq!(ldexp(-1.25, 0), -1.25);
    }

    #[test]
    fn coefficient_scaling_keeps_average_finite_and_total_infinite() {
        // The product, the sum of magnitudes, and a single coefficient
        // respectively exceed binary64 at these three sizes.
        for (n, expected) in [
            (1473, 5.567935488305321),
            (1476, 5.569403103260331),
            (1500, 5.581037752840027),
        ] {
            assert_eq!(information_content(path_graph(n), false), f64::INFINITY);
            assert_eq!(information_content(path_graph(n), true), expected);
        }
    }

    #[test]
    fn diagonal_loop_is_included_in_the_characteristic_polynomial() {
        let mut looped = graph(1, &[]);
        looped.loops[0] = true;
        assert_eq!(
            characteristic_polynomial_magnitudes(&looped)
                .unwrap()
                .scaled,
            [1.0, 1.0]
        );
        let mut looped_again = graph(1, &[]);
        looped_again.loops[0] = true;
        assert_eq!(information_content(looped_again, false), 2.0);
        let mut looped_average = graph(1, &[]);
        looped_average.loops[0] = true;
        assert_eq!(information_content(looped_average, true), 1.0);
    }

    fn binomial(n: usize, k: usize) -> u128 {
        let k = k.min(n - k);
        (0..k).fold(1_u128, |value, i| value * (n - i) as u128 / (i + 1) as u128)
    }
}
