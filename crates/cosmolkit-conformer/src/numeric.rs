//! Complete source distance sampling, eigen initial coordinates and random coordinates.
use crate::{ConformerError, bounds::BoundsMatrix};
use cosmolkit_core::{RdkitRandomEngine, RdkitRandomGenerator, with_rdkit_random_generator};
const EIGVAL_TOL: f64 = 0.001;
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct SymmMatrix<T = f64> {
    data: Vec<T>,
    n: usize,
}

impl<T: Copy + Default> SymmMatrix<T> {
    pub(crate) fn new(n: usize) -> Self {
        Self {
            data: vec![T::default(); n * (n + 1) / 2],
            n,
        }
    }

    pub(crate) fn with_value(n: usize, value: T) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR RDNumeric::SymmMatrix::SymmMatrix value overload (SymmMatrix.h:40-48)
        // RDKit✔️✔️:   SymmMatrix(unsigned int N, TYPE val)
        // RDKit✔️✔️:       : d_size(N), d_dataSize(N * (N + 1) / 2) {
        // RDKit✔️✔️:     TYPE *data = new TYPE[d_dataSize];
        // RDKit✔️✔️:     unsigned int i;
        // RDKit✔️✔️:     for (i = 0; i < d_dataSize; i++) {
        // RDKit✔️✔️:       data[i] = val;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     d_data.reset(data);
        // RDKit✔️✔️:   }
        // END RDKIT CPP CONSTRUCTOR RDNumeric::SymmMatrix::SymmMatrix value overload
        Self {
            data: vec![value; n * (n + 1) / 2],
            n,
        }
    }

    pub(crate) fn num_rows(&self) -> usize {
        self.n
    }

    pub(crate) fn get_data_size(&self) -> usize {
        self.data.len()
    }

    pub(crate) fn get_data(&self) -> &[T] {
        &self.data
    }

    pub(crate) fn get_val(&self, i: usize, j: usize) -> T {
        self.data[Self::data_index(i, j)]
    }

    pub(crate) fn set_val(&mut self, i: usize, j: usize, value: T) {
        let idx = Self::data_index(i, j);
        self.data[idx] = value;
    }

    pub(crate) fn data_index(i: usize, j: usize) -> usize {
        let (row, col) = if i >= j { (i, j) } else { (j, i) };
        row * (row + 1) / 2 + col
    }
}

#[derive(Debug, Clone, PartialEq)]
struct DoubleMatrix {
    data: Vec<f64>,
    n_rows: usize,
    n_cols: usize,
}

impl DoubleMatrix {
    fn new(n_rows: usize, n_cols: usize) -> Self {
        // BEGIN RDKIT CPP CONSTRUCTOR RDNumeric::Matrix::Matrix size overload (Matrix.h:33-38)
        // RDKit✔️✔️:   Matrix(unsigned int nRows, unsigned int nCols)
        // RDKit✔️✔️:       : d_nRows(nRows), d_nCols(nCols), d_dataSize(nRows * nCols) {
        // RDKit✔️✔️:     TYPE *data = new TYPE[d_dataSize];
        // RDKit✔️✔️:     memset(static_cast<void *>(data), 0, d_dataSize * sizeof(TYPE));
        // RDKit✔️✔️:     d_data.reset(data);
        // RDKit✔️✔️:   }
        // END RDKIT CPP CONSTRUCTOR RDNumeric::Matrix::Matrix size overload
        Self {
            data: vec![0.0; n_rows * n_cols],
            n_rows,
            n_cols,
        }
    }

    fn num_rows(&self) -> usize {
        self.n_rows
    }

    fn num_cols(&self) -> usize {
        self.n_cols
    }

    fn get_val(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.n_cols + j]
    }
}

pub(crate) trait RdkitDoubleRng {
    fn next_unit_f64(&mut self) -> f64;
}

#[cfg(not(target_arch = "wasm32"))]
unsafe extern "C" {
    fn clock() -> libc::clock_t;
}

#[cfg(not(target_arch = "wasm32"))]
pub(crate) fn rdkit_clock_seed() -> Result<i32, ConformerError> {
    // BEGIN RDKIT CPP CLOCK SEED SOURCE Vector.h:279 and EigenSolvers/PowerEigenSolver.cpp:47
    // RDKit✔️✔️:       generator.seed(clock() + 1);
    // RDKit✔️✔️:     seed = clock();
    // END RDKIT CPP CLOCK SEED SOURCE
    // Native source callers apply their own +1 and signed narrowing.
    Ok(unsafe { clock() as i32 })
}

#[cfg(target_arch = "wasm32")]
pub(crate) fn rdkit_clock_seed() -> Result<i32, ConformerError> {
    // BEGIN RDKIT CPP CLOCK SEED SOURCE Vector.h:279 and EigenSolvers/PowerEigenSolver.cpp:47
    // RDKit❌❌:       generator.seed(clock() + 1);
    // RDKit❌❌:     seed = clock();
    // END RDKIT CPP CLOCK SEED SOURCE
    // wasm32 has no modeled C process clock; the independent implicit-seed
    // capability remains a structured unsupported error, with no substitute.
    Err(ConformerError::WasmImplicitClockSeedUnsupported)
}

impl RdkitDoubleRng for RdkitRandomEngine {
    fn next_unit_f64(&mut self) -> f64 {
        RdkitRandomEngine::next_unit_f64(self)
    }
}
impl RdkitDoubleRng for RdkitRandomGenerator<'_> {
    fn next_unit_f64(&mut self) -> f64 {
        RdkitRandomGenerator::next_unit_f64(self)
    }
}
pub(crate) fn rdkit_embedder_multiplication_overflows(a: i32, b: i32) -> bool {
    // BEGIN RDKIT CPP TEMPLATE FUNCTION DGeomHelpers::detail::multiplication_overflows_ (Embedder.cpp:1346-1353)
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: bool multiplication_overflows_(T a, T b) {
    // RDKit✔️✔️:   // a * b > c if and only if a > c / b
    // RDKit✔️✔️:   if (a == 0 || b == 0) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return a > std::numeric_limits<T>::max() / b;
    // RDKit✔️✔️: }
    // END RDKIT CPP TEMPLATE FUNCTION DGeomHelpers::detail::multiplication_overflows_
    if a == 0 || b == 0 {
        return false;
    }
    a > i32::MAX / b
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum EmbedSeedError {
    #[error("random seed must either be positive, zero, or negative one")]
    InvalidRandomSeed { seed: i32 },
    #[error("Something went wrong calculating a new seed")]
    ResultInvariant { seed: i32, index: usize, value: i32 },
}
pub(crate) fn rdkit_embedder_conformer_seed(
    random_seed: i32,
    conformer_index: usize,
    enable_sequential_random_seeds: bool,
) -> Result<i32, EmbedSeedError> {
    // BEGIN RDKIT CPP MACRO release rdcast (RDGeneral/Invariant.h:183-193)
    // RDKit❗✔️: #ifdef RDDEBUG
    // RDKit❗✔️: // use rdcast to convert between types
    // RDKit❗✔️: //  when RDDEBUG is defined, this checks for
    // RDKit❗✔️: //  validity (overflow, etc)
    // RDKit❗✔️: //  when RDDEBUG is off, the cast is a no-cost
    // RDKit❗✔️: //   static_cast
    // RDKit❗✔️: #define rdcast boost::numeric_cast
    // RDKit❗✔️: #else
    // RDKit❗✔️: #define rdcast static_cast
    // RDKit❗✔️: #endif
    // END RDKIT CPP MACRO release rdcast

    // BEGIN RDKIT CPP BLOCK DGeomHelpers::detail::embedHelper_ per-conformer seed policy (Embedder.cpp:1393-1426)
    // RDKit❗✔️:     CHECK_INVARIANT(
    // RDKit❗✔️:         params->randomSeed >= -1,
    // RDKit❗✔️:         "random seed must either be positive, zero, or negative one");
    // RDKit❗✔️:     int new_seed = params->randomSeed;
    // RDKit❗✔️:     if (new_seed > -1) {
    // RDKit❗✔️:       if (params->enableSequentialRandomSeeds) {
    // RDKit❗✔️:         new_seed += ci + 1;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         if (!multiplication_overflows_(rdcast<int>(ci + 1),
    // RDKit❗✔️:                                        params->randomSeed)) {
    // RDKit❗✔️:           // old method of computing a new seed
    // RDKit❗✔️:           new_seed = (ci + 1) * params->randomSeed;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           // If the above simple multiplication will overflow, use a
    // RDKit❗✔️:           // cheap and easy way to hash the conformer index and seed
    // RDKit❗✔️:           // together: for N'ary numerical system, where N is the
    // RDKit❗✔️:           // maximum possible value of the pair of numbers. The
    // RDKit❗✔️:           // following will generate unique integers:
    // RDKit❗✔️:           // hash(a, b) = a + b * N
    // RDKit❗✔️:           auto big_seed = rdcast<size_t>(params->randomSeed);
    // RDKit❗✔️:           size_t max_val = std::max(ci + 1, big_seed);
    // RDKit❗✔️:           size_t big_num = big_seed + max_val * (ci + 1);
    // RDKit❗✔️:           // only grab the first 31 bits xor'd with the next 31 bits to
    // RDKit❗✔️:           // make sure its positive, careful, the 'ULL' is important
    // RDKit❗✔️:           // here, 0x7fffffff is the 'int' type because of C default
    // RDKit❗✔️:           // number semantics and that we definitely don't want!
    // RDKit❗✔️:           const size_t positive_int_mask = 0x7fffffffULL;
    // RDKit❗✔️:           size_t folded_num =
    // RDKit❗✔️:               (big_num & positive_int_mask) ^ (big_num >> 31ULL);
    // RDKit❗✔️:           new_seed = rdcast<int>(folded_num & positive_int_mask);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     CHECK_INVARIANT(new_seed >= -1,
    // RDKit❗✔️:                     "Something went wrong calculating a new seed");
    // END RDKIT CPP BLOCK DGeomHelpers::detail::embedHelper_ per-conformer seed policy
    // The pinned release rdcast is static_cast (RDDEBUG is a distinct
    // upstream diagnostic build). Compound int += size_t and size_t products
    // therefore wrap at pointer width and narrow to the source signed int.
    if random_seed < -1 {
        return Err(EmbedSeedError::InvalidRandomSeed { seed: random_seed });
    }
    let mut new_seed = random_seed;
    if new_seed > -1 {
        let index = conformer_index.wrapping_add(1);
        if enable_sequential_random_seeds {
            new_seed = (new_seed as usize).wrapping_add(index) as i32;
        } else if !rdkit_embedder_multiplication_overflows(index as i32, random_seed) {
            new_seed = index.wrapping_mul(random_seed as usize) as i32;
        } else {
            let big_seed = random_seed as usize;
            let max_val = index.max(big_seed);
            let big_num = big_seed.wrapping_add(max_val.wrapping_mul(index));
            const POSITIVE_INT_MASK: usize = 0x7fff_ffff;
            let folded = (big_num & POSITIVE_INT_MASK) ^ (big_num >> 31);
            new_seed = (folded & POSITIVE_INT_MASK) as i32;
        }
    }
    if new_seed < -1 {
        return Err(EmbedSeedError::ResultInvariant {
            seed: random_seed,
            index: conformer_index,
            value: new_seed,
        });
    }
    Ok(new_seed)
}

pub(crate) fn pick_random_dist_mat(
    mmat: &BoundsMatrix,
    dist_mat: &mut SymmMatrix,
    seed: i32,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::pickRandomDistMat seed overload (DistGeomUtils.cpp:35-41)
    // RDKit✔️✔️: double pickRandomDistMat(const BoundsMatrix &mmat,
    // RDKit✔️✔️:                          RDNumeric::SymmMatrix<double> &distMat, int seed) {
    // RDKit✔️✔️:   if (seed > 0) {
    // RDKit✔️✔️:     RDKit::getRandomGenerator(seed);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return pickRandomDistMat(mmat, distMat, RDKit::getDoubleRandomSource());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::pickRandomDistMat seed overload
    // One borrow of the existing source process stream for the full operation.
    with_rdkit_random_generator(seed, |rng| {
        pick_random_dist_mat_with_rng(mmat, dist_mat, rng)
    })
}

pub(crate) fn pick_random_dist_mat_with_rng<R: RdkitDoubleRng>(
    mmat: &BoundsMatrix,
    dist_mat: &mut SymmMatrix,
    rng: &mut R,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION DistGeom::pickRandomDistMat RNG overload (DistGeomUtils.cpp:43-68)
    // RDKit✔️✔️: double pickRandomDistMat(const BoundsMatrix &mmat,
    // RDKit✔️✔️:                          RDNumeric::SymmMatrix<double> &distMat,
    // RDKit✔️✔️:                          RDKit::double_source_type &rng) {
    // RDKit✔️✔️:   // make sure the sizes match up
    // RDKit✔️✔️:   unsigned int npt = mmat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(npt == distMat.numRows(), "Size mismatch");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double largestVal = -1.0;
    // RDKit✔️✔️:   double *ddata = distMat.getData();
    // RDKit✔️✔️:   for (unsigned int i = 1; i < npt; i++) {
    // RDKit✔️✔️:     unsigned int id = i * (i + 1) / 2;
    // RDKit✔️✔️:     for (unsigned int j = 0; j < i; j++) {
    // RDKit✔️✔️:       double ub = mmat.getUpperBound(i, j);
    // RDKit✔️✔️:       double lb = mmat.getLowerBound(i, j);
    // RDKit✔️✔️:       CHECK_INVARIANT(ub >= lb, "");
    // RDKit✔️✔️:       double rval = rng();
    // RDKit✔️✔️:       // std::cerr<<i<<"-"<<j<<": "<<rval<<std::endl;
    // RDKit✔️✔️:       double d = lb + (rval) * (ub - lb);
    // RDKit✔️✔️:       ddata[id + j] = d;
    // RDKit✔️✔️:       if (d > largestVal) {
    // RDKit✔️✔️:         largestVal = d;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return largestVal;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::pickRandomDistMat RNG overload
    let npt = mmat.dimension();
    assert_eq!(npt, dist_mat.num_rows(), "Size mismatch");

    let mut largest_val = -1.0;
    for i in 1..npt {
        let id = i * (i + 1) / 2;
        for j in 0..i {
            let ub = mmat.get_upper(i, j).expect("source upper-bound indices");
            let lb = mmat.get_lower(i, j).expect("source lower-bound indices");
            assert!(ub >= lb);
            let rval = rng.next_unit_f64();
            let d = lb + rval * (ub - lb);
            dist_mat.data[id + j] = d;
            if d > largest_val {
                largest_val = d;
            }
        }
    }

    largest_val
}

pub(crate) fn rdkit_symm_matrix_vector_multiply(a: &SymmMatrix, x: &[f64], y: &mut [f64]) {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::multiply SymmMatrix-Vector overload (SymmMatrix.h:313-335)
    // RDKit✔️✔️: Vector<TYPE> &multiply(const SymmMatrix<TYPE> &A, const Vector<TYPE> &x,
    // RDKit✔️✔️:                        Vector<TYPE> &y) {
    // RDKit✔️✔️:   unsigned int aSize = A.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(aSize == x.size(), "Size mismatch during multiplication");
    // RDKit✔️✔️:   CHECK_INVARIANT(aSize == y.size(), "Size mismatch during multiplication");
    // RDKit✔️✔️:   const TYPE *xData = x.getData();
    // RDKit✔️✔️:   const TYPE *aData = A.getData();
    // RDKit✔️✔️:   TYPE *yData = y.getData();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < aSize; i++) {
    // RDKit✔️✔️:     yData[i] = (TYPE)(0.0);
    // RDKit✔️✔️:     unsigned int idA = i * (i + 1) / 2;
    // RDKit✔️✔️:     for (unsigned int j = 0; j < i + 1; j++) {
    // RDKit✔️✔️:       // idA = i*(i+1)/2 + j;
    // RDKit✔️✔️:       yData[i] += (aData[idA] * xData[j]);
    // RDKit✔️✔️:       idA++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     idA--;
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < aSize; j++) {
    // RDKit✔️✔️:       // idA = j*(j+1)/2 + i;
    // RDKit✔️✔️:       idA += j;
    // RDKit✔️✔️:       yData[i] += (aData[idA] * xData[j]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDNumeric::multiply SymmMatrix-Vector overload
    let a_size = a.num_rows();
    assert_eq!(a_size, x.len(), "Size mismatch during multiplication");
    assert_eq!(a_size, y.len(), "Size mismatch during multiplication");
    for i in 0..a_size {
        y[i] = 0.0;
        let mut id_a = i * (i + 1) / 2;
        for xj in x.iter().take(i + 1) {
            y[i] += a.data[id_a] * *xj;
            id_a += 1;
        }
        id_a -= 1;
        for (j, xj) in x.iter().enumerate().take(a_size).skip(i + 1) {
            id_a += j;
            y[i] += a.data[id_a] * *xj;
        }
    }
}

pub(crate) fn rdkit_vector_largest_abs_val_id(data: &[f64]) -> usize {
    // BEGIN RDKIT CPP METHOD RDNumeric::Vector::largestAbsValId (Vector.h:202-213)
    // RDKit✔️✔️:   constexpr unsigned int largestAbsValId() const {
    // RDKit✔️✔️:     TYPE res = (TYPE)(-1.0);
    // RDKit✔️✔️:     unsigned int i, id = d_size;
    // RDKit✔️✔️:     TYPE *data = d_data.get();
    // RDKit✔️✔️:     for (i = 0; i < d_size; i++) {
    // RDKit✔️✔️:       if (fabs(data[i]) > res) {
    // RDKit✔️✔️:         res = fabs(data[i]);
    // RDKit✔️✔️:         id = i;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return id;
    // RDKit✔️✔️:   }
    // END RDKIT CPP METHOD RDNumeric::Vector::largestAbsValId
    let mut res = -1.0;
    let mut id = data.len();
    for (i, value) in data.iter().enumerate() {
        if value.abs() > res {
            res = value.abs();
            id = i;
        }
    }
    id
}

pub(crate) fn rdkit_vector_normalize(data: &mut [f64]) -> Result<(), ConformerError> {
    // RDKit✔️✔️: static constexpr double zero_tolerance = 1.e-16;
    // BEGIN RDKIT CPP METHOD RDNumeric::Vector::normalize (Vector.h:258-264)
    // RDKit✔️✔️:   constexpr void normalize() {
    // RDKit✔️✔️:     TYPE val = this->normL2();
    // RDKit✔️✔️:     if (val < zero_tolerance) {
    // RDKit✔️✔️:       throw std::runtime_error("Cannot normalize a zero length vector");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     (*this) /= val;
    // RDKit✔️✔️:   }
    // END RDKIT CPP METHOD RDNumeric::Vector::normalize
    let val = data.iter().map(|v| v * v).sum::<f64>().sqrt();
    if val < 1.0e-16 {
        return Err(ConformerError::CannotNormalizeZeroLengthVector);
    }
    for item in data {
        *item /= val;
    }
    Ok(())
}

pub(crate) fn rdkit_vector_set_to_random(
    size: usize,
    seed: i32,
) -> Result<Vec<f64>, ConformerError> {
    // BEGIN RDKIT CPP METHOD RDNumeric::Vector::setToRandom (Vector.h:267-288)
    // RDKit✔️✔️:   void setToRandom(unsigned int seed = 0) {
    // RDKit✔️✔️:     // we want to get our own RNG here instead of using the global
    // RDKit✔️✔️:     // one.  This is related to Issue285.
    // RDKit✔️✔️:     RDKit::rng_type generator(42u);
    // RDKit✔️✔️:     RDKit::uniform_double dist(0, 1.0);
    // RDKit✔️✔️:     RDKit::double_source_type randSource(generator, dist);
    // RDKit✔️✔️:     if (seed > 0) {
    // RDKit✔️✔️:       generator.seed(seed);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // we can't initialize using only clock(), because it's possible
    // RDKit✔️✔️:       // that we'll get here fast enough that clock() will return 0
    // RDKit✔️✔️:       // and generator.seed(0) is an error:
    // RDKit✔️✔️:       generator.seed(clock() + 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     unsigned int i;
    // RDKit✔️✔️:     TYPE *data = d_data.get();
    // RDKit✔️✔️:     for (i = 0; i < d_size; i++) {
    // RDKit✔️✔️:       data[i] = randSource();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     this->normalize();
    // RDKit✔️✔️:   }
    // END RDKIT CPP METHOD RDNumeric::Vector::setToRandom
    let effective_seed = if seed > 0 {
        seed
    } else {
        rdkit_clock_seed()?.wrapping_add(1)
    };
    let mut rng = RdkitRandomEngine::from_seed(effective_seed as u32);
    let mut data = (0..size).map(|_| rng.next_unit_f64()).collect::<Vec<_>>();
    rdkit_vector_normalize(&mut data)?;
    Ok(data)
}

fn power_eigen_solver(
    num_eig: usize,
    mat: &mut SymmMatrix,
    eigen_values: &mut [f64],
    mut eigen_vectors: Option<&mut DoubleMatrix>,
    seed: i32,
) -> Result<bool, ConformerError> {
    // BEGIN RDKIT CPP FUNCTION RDNumeric::EigenSolvers::powerEigenSolver (PowerEigenSolver.cpp:20-100)
    // RDKit✔️✔️: bool powerEigenSolver(unsigned int numEig, DoubleSymmMatrix &mat,
    // RDKit✔️✔️:                       DoubleVector &eigenValues, DoubleMatrix *eigenVectors,
    // RDKit✔️✔️:                       int seed) {
    // RDKit✔️✔️:   const unsigned int MAX_ITERATIONS = 1000;
    // RDKit✔️✔️:   const double TOLERANCE = 0.001;
    // RDKit✔️✔️:   const double HUGE_EIGVAL = 1.0e10;
    // RDKit✔️✔️:   const double TINY_EIGVAL = 1.0e-10;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // first check all the sizes
    // RDKit✔️✔️:   unsigned int N = mat.numRows();
    // RDKit✔️✔️:   CHECK_INVARIANT(eigenValues.size() >= numEig, "");
    // RDKit✔️✔️:   CHECK_INVARIANT(numEig <= N, "");
    // RDKit✔️✔️:   if (eigenVectors) {
    // RDKit✔️✔️:     unsigned int evRows, evCols;
    // RDKit✔️✔️:     evRows = eigenVectors->numRows();
    // RDKit✔️✔️:     evCols = eigenVectors->numCols();
    // RDKit✔️✔️:     CHECK_INVARIANT(evCols >= N, "");
    // RDKit✔️✔️:     CHECK_INVARIANT(evRows >= numEig, "");
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION RDNumeric::EigenSolvers::powerEigenSolver
    const MAX_ITERATIONS: usize = 1000;
    const TOLERANCE: f64 = 0.001;
    const HUGE_EIGVAL: f64 = 1.0e10;
    const TINY_EIGVAL: f64 = 1.0e-10;

    let n = mat.num_rows();
    assert!(eigen_values.len() >= num_eig);
    assert!(num_eig <= n);
    if let Some(eig_vecs) = eigen_vectors.as_ref() {
        assert!(eig_vecs.num_cols() >= n);
        assert!(eig_vecs.num_rows() >= num_eig);
    }

    // RDKit✔️✔️:   unsigned int ei;
    // RDKit✔️✔️:   double eigVal, prevVal;
    // RDKit✔️✔️:   bool converged = false;
    // RDKit✔️✔️:   unsigned int i, j, id, iter, evalId;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   DoubleVector v(N), z(N);
    // RDKit✔️✔️:   if (seed <= 0) {
    // RDKit✔️✔️:     seed = clock();
    // RDKit✔️✔️:   }
    let mut effective_seed = if seed > 0 { seed } else { rdkit_clock_seed()? };
    let mut converged = false;
    let mut z = vec![0.0; n];
    for (ei, eig_value_slot) in eigen_values.iter_mut().enumerate().take(num_eig) {
        // RDKit✔️✔️:   for (ei = 0; ei < numEig; ei++) {
        // RDKit✔️✔️:     eigVal = -HUGE_EIGVAL;
        // RDKit✔️✔️:     seed += ei;
        // RDKit✔️✔️:     v.setToRandom(seed);
        let mut eig_val = -HUGE_EIGVAL;
        effective_seed += ei as i32;
        let mut v = rdkit_vector_set_to_random(n, effective_seed)?;

        converged = false;
        for _iter in 0..MAX_ITERATIONS {
            // RDKit✔️✔️:     for (iter = 0; iter < MAX_ITERATIONS; iter++) {
            // RDKit✔️✔️:       // z = mat*v
            // RDKit✔️✔️:       multiply(mat, v, z);
            // RDKit✔️✔️:       prevVal = eigVal;
            // RDKit✔️✔️:       evalId = z.largestAbsValId();
            // RDKit✔️✔️:       eigVal = z.getVal(evalId);
            rdkit_symm_matrix_vector_multiply(mat, &v, &mut z);
            let prev_val = eig_val;
            let eval_id = rdkit_vector_largest_abs_val_id(&z);
            eig_val = z[eval_id];

            // RDKit✔️✔️:       if (fabs(eigVal) < TINY_EIGVAL) {
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            if eig_val.abs() < TINY_EIGVAL {
                break;
            }

            // RDKit✔️✔️:       // compute the next estimate for the eigen vector
            // RDKit✔️✔️:       v.assign(z);
            // RDKit✔️✔️:       v /= eigVal;
            // RDKit✔️✔️:       if (fabs(eigVal - prevVal) < TOLERANCE) {
            // RDKit✔️✔️:         converged = true;
            // RDKit✔️✔️:         break;
            // RDKit✔️✔️:       }
            v.copy_from_slice(&z);
            for item in &mut v {
                *item /= eig_val;
            }
            if (eig_val - prev_val).abs() < TOLERANCE {
                converged = true;
                break;
            }
        }
        // RDKit✔️✔️:     if (!converged) {
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        if !converged {
            break;
        }
        // RDKit✔️✔️:     v.normalize();
        rdkit_vector_normalize(&mut v)?;

        // RDKit✔️✔️:     // save this is a eigen vector and value
        // RDKit✔️✔️:     // directly access the data instead of setVal so that we save time
        // RDKit✔️✔️:     double *vdata = v.getData();
        // RDKit✔️✔️:     if (eigenVectors) {
        // RDKit✔️✔️:       id = ei * eigenVectors->numCols();
        // RDKit✔️✔️:       double *eigVecData = eigenVectors->getData();
        // RDKit✔️✔️:       for (i = 0; i < N; i++) {
        // RDKit✔️✔️:         eigVecData[id + i] = vdata[i];
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        if let Some(eig_vecs) = eigen_vectors.as_deref_mut() {
            let id = ei * eig_vecs.num_cols();
            eig_vecs.data[id..id + n].copy_from_slice(&v);
        }
        // RDKit✔️✔️:     eigenValues[ei] = eigVal;
        *eig_value_slot = eig_val;

        // RDKit✔️✔️:     // now remove this eigen vector space out of the matrix
        // RDKit✔️✔️:     double *matData = mat.getData();
        // RDKit✔️✔️:     for (i = 0; i < N; i++) {
        // RDKit✔️✔️:       id = i * (i + 1) / 2;
        // RDKit✔️✔️:       for (j = 0; j < i + 1; j++) {
        // RDKit✔️✔️:         matData[id + j] -= (eigVal * vdata[i] * vdata[j]);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for i in 0..n {
            let id = i * (i + 1) / 2;
            for j in 0..=i {
                mat.data[id + j] -= eig_val * v[i] * v[j];
            }
        }
    }
    // RDKit✔️✔️:   return converged;
    // RDKit✔️✔️: }
    Ok(converged)
}

pub(crate) fn compute_initial_coords(
    dist_mat: &SymmMatrix,
    positions: &mut [Vec<f64>],
    rand_neg_eig: bool,
    num_zero_fail: usize,
    seed: i32,
) -> Result<bool, ConformerError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeInitialCoords seed overload (DistGeomUtils.cpp:70-80)
    // RDKit✔️✔️: bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,
    // RDKit✔️✔️:                           RDGeom::PointPtrVect &positions, bool randNegEig,
    // RDKit✔️✔️:                           unsigned int numZeroFail, int seed) {
    // RDKit✔️✔️:   if (seed > 0) {
    // RDKit✔️✔️:     RDKit::getRandomGenerator(seed);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return computeInitialCoords(distMat, positions,
    // RDKit✔️✔️:                               RDKit::getDoubleRandomSource(), randNegEig,
    // RDKit✔️✔️:                               numZeroFail);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::computeInitialCoords seed overload
    // One borrow of the existing source process stream for the full operation.
    with_rdkit_random_generator(seed, |rng| {
        compute_initial_coords_with_rng(dist_mat, positions, rng, rand_neg_eig, num_zero_fail)
    })
}

pub(crate) fn compute_initial_coords_with_rng<R: RdkitDoubleRng>(
    dist_mat: &SymmMatrix,
    positions: &mut [Vec<f64>],
    rng: &mut R,
    rand_neg_eig: bool,
    num_zero_fail: usize,
) -> Result<bool, ConformerError> {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeInitialCoords RNG overload (DistGeomUtils.cpp:81-164)
    // RDKit✔️✔️: bool computeInitialCoords(const RDNumeric::SymmMatrix<double> &distMat,
    // RDKit✔️✔️:                           RDGeom::PointPtrVect &positions,
    // RDKit✔️✔️:                           RDKit::double_source_type &rng, bool randNegEig,
    // RDKit✔️✔️:                           unsigned int numZeroFail) {
    // RDKit✔️✔️:   unsigned int N = distMat.numRows();
    // RDKit✔️✔️:   unsigned int nPt = positions.size();
    // RDKit✔️✔️:   CHECK_INVARIANT(nPt == N, "Size mismatch");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int dim = positions.front()->dimension();
    // END RDKIT CPP FUNCTION DistGeom::computeInitialCoords RNG overload
    let n = dist_mat.num_rows();
    assert_eq!(positions.len(), n, "Size mismatch");
    assert!(!positions.is_empty());
    let dim = positions[0].len();

    // RDKit✔️✔️:   const double *data = distMat.getData();
    // RDKit✔️✔️:   RDNumeric::SymmMatrix<double> sqMat(N), T(N, 0.0);
    // RDKit✔️✔️:   RDNumeric::DoubleMatrix eigVecs(dim, N);
    // RDKit✔️✔️:   RDNumeric::DoubleVector eigVals(dim);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   double *sqDat = sqMat.getData();
    let data = dist_mat.get_data();
    let mut sq_mat = SymmMatrix::new(n);
    let mut t = SymmMatrix::with_value(n, 0.0);
    let mut eig_vecs = DoubleMatrix::new(dim, n);
    let mut eig_vals = vec![0.0; dim];

    // RDKit✔️✔️:   unsigned int dSize = distMat.getDataSize();
    // RDKit✔️✔️:   double sumSqD2 = 0.0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dSize; i++) {
    // RDKit✔️✔️:     sqDat[i] = data[i] * data[i];
    // RDKit✔️✔️:     sumSqD2 += sqDat[i];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   sumSqD2 /= (N * N);
    let d_size = dist_mat.get_data_size();
    let mut sum_sq_d2 = 0.0;
    for (i, value) in data.iter().enumerate().take(d_size) {
        sq_mat.data[i] = value * value;
        sum_sq_d2 += sq_mat.data[i];
    }
    sum_sq_d2 /= (n * n) as f64;

    // RDKit✔️✔️:   RDNumeric::DoubleVector sqD0i(N, 0.0);
    // RDKit✔️✔️:   double *sqD0iData = sqD0i.getData();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:     for (unsigned int j = 0; j < N; j++) {
    // RDKit✔️✔️:       sqD0iData[i] += sqMat.getVal(i, j);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     sqD0iData[i] /= N;
    // RDKit✔️✔️:     sqD0iData[i] -= sumSqD2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     if ((sqD0iData[i] < EIGVAL_TOL) && (N > 3)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut sq_d0i = vec![0.0; n];
    for (i, sq_d0i_value) in sq_d0i.iter_mut().enumerate().take(n) {
        for j in 0..n {
            *sq_d0i_value += sq_mat.get_val(i, j);
        }
        *sq_d0i_value /= n as f64;
        *sq_d0i_value -= sum_sq_d2;
        if *sq_d0i_value < EIGVAL_TOL && n > 3 {
            return Ok(false);
        }
    }

    // RDKit✔️✔️:   for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:     for (unsigned int j = 0; j <= i; j++) {
    // RDKit✔️✔️:       double val = 0.5 * (sqD0iData[i] + sqD0iData[j] - sqMat.getVal(i, j));
    // RDKit✔️✔️:       T.setVal(i, j, val);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int nEigs = (dim < N) ? dim : N;
    // RDKit✔️✔️:   RDNumeric::EigenSolvers::powerEigenSolver(nEigs, T, eigVals, eigVecs,
    // RDKit✔️✔️:                                             (int)(sumSqD2 * N));
    for i in 0..n {
        for j in 0..=i {
            let val = 0.5 * (sq_d0i[i] + sq_d0i[j] - sq_mat.get_val(i, j));
            t.set_val(i, j, val);
        }
    }
    let n_eigs = dim.min(n);
    power_eigen_solver(
        n_eigs,
        &mut t,
        &mut eig_vals,
        Some(&mut eig_vecs),
        (sum_sq_d2 * n as f64) as i32,
    )?;

    // RDKit✔️✔️:   double *eigData = eigVals.getData();
    // RDKit✔️✔️:   bool foundNeg = false;
    // RDKit✔️✔️:   unsigned int zeroEigs = 0;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < dim; i++) {
    // RDKit✔️✔️:     if (eigData[i] > EIGVAL_TOL) {
    // RDKit✔️✔️:       eigData[i] = sqrt(eigData[i]);
    // RDKit✔️✔️:     } else if (fabs(eigData[i]) < EIGVAL_TOL) {
    // RDKit✔️✔️:       eigData[i] = 0.0;
    // RDKit✔️✔️:       zeroEigs++;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       foundNeg = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut found_neg = false;
    let mut zero_eigs = 0;
    for eig_val in eig_vals.iter_mut().take(dim) {
        if *eig_val > EIGVAL_TOL {
            *eig_val = eig_val.sqrt();
        } else if eig_val.abs() < EIGVAL_TOL {
            *eig_val = 0.0;
            zero_eigs += 1;
        } else {
            found_neg = true;
        }
    }

    // RDKit✔️✔️:   if ((foundNeg) && (!randNegEig)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if ((zeroEigs >= numZeroFail) && (N > 3)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    if found_neg && !rand_neg_eig {
        return Ok(false);
    }
    if zero_eigs >= num_zero_fail && n > 3 {
        return Ok(false);
    }

    // RDKit✔️✔️:   for (unsigned int i = 0; i < N; i++) {
    // RDKit✔️✔️:     RDGeom::Point *pt = positions[i];
    // RDKit✔️✔️:     for (unsigned int j = 0; j < dim; ++j) {
    // RDKit✔️✔️:       if (eigData[j] >= 0.0) {
    // RDKit✔️✔️:         (*pt)[j] = eigData[j] * eigVecs.getVal(j, i);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         // std::cerr<<"!!! "<<i<<"-"<<j<<": "<<eigData[j]<<std::endl;
    // RDKit✔️✔️:         (*pt)[j] = 1.0 - 2.0 * rng();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    for (i, point) in positions.iter_mut().enumerate().take(n) {
        for j in 0..dim {
            if eig_vals[j] >= 0.0 {
                point[j] = eig_vals[j] * eig_vecs.get_val(j, i);
            } else {
                point[j] = 1.0 - 2.0 * rng.next_unit_f64();
            }
        }
    }
    Ok(true)
}

pub(crate) fn compute_random_coords(positions: &mut [Vec<f64>], box_size: f64, seed: i32) -> bool {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeRandomCoords seed overload (DistGeomUtils.cpp:166-172)
    // RDKit✔️✔️: bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
    // RDKit✔️✔️:                          int seed) {
    // RDKit✔️✔️:   if (seed > 0) {
    // RDKit✔️✔️:     RDKit::getRandomGenerator(seed);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return computeRandomCoords(positions, boxSize,
    // RDKit✔️✔️:                              RDKit::getDoubleRandomSource());
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::computeRandomCoords seed overload
    // One borrow of the existing source process stream for the full operation.
    with_rdkit_random_generator(seed, |rng| {
        compute_random_coords_with_rng(positions, box_size, rng)
    })
}

pub(crate) fn compute_random_coords_with_rng<R: RdkitDoubleRng>(
    positions: &mut [Vec<f64>],
    box_size: f64,
    rng: &mut R,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION DistGeom::computeRandomCoords RNG overload (DistGeomUtils.cpp:173-183)
    // RDKit✔️✔️: bool computeRandomCoords(RDGeom::PointPtrVect &positions, double boxSize,
    // RDKit✔️✔️:                          RDKit::double_source_type &rng) {
    // RDKit✔️✔️:   CHECK_INVARIANT(boxSize > 0.0, "bad boxSize");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (auto pt : positions) {
    // RDKit✔️✔️:     for (unsigned int i = 0; i < pt->dimension(); ++i) {
    // RDKit✔️✔️:       (*pt)[i] = boxSize * (rng() - 0.5);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION DistGeom::computeRandomCoords RNG overload
    assert!(box_size > 0.0, "bad boxSize");
    for point in positions {
        for coord in point {
            *coord = box_size * (rng.next_unit_f64() - 0.5);
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;
    // Proposed regressions for pinned Vector.h:24,258-264. Original parity
    // conditions are unchanged; independent review is recorded separately.
    #[test]
    fn vector_normalize_source_below_threshold_preserves_input_on_error() {
        let mut data = [0.5e-16, -0.0];
        let before = data.map(f64::to_bits);
        let error = rdkit_vector_normalize(&mut data).unwrap_err();
        assert_eq!(error, ConformerError::CannotNormalizeZeroLengthVector);
        assert_eq!(error.to_string(), "Cannot normalize a zero length vector");
        assert_eq!(data.map(f64::to_bits), before);
    }

    #[test]
    fn vector_normalize_source_equal_threshold() {
        let mut data = [1.0e-16];
        rdkit_vector_normalize(&mut data).unwrap();
        assert_eq!(data, [1.0]);
    }

    #[test]
    fn vector_normalize_source_between_threshold_and_epsilon() {
        let mut data = [1.5e-16];
        rdkit_vector_normalize(&mut data).unwrap();
        assert_eq!(data, [1.0]);
    }

    #[test]
    fn vector_normalize_source_epsilon() {
        let mut data = [f64::EPSILON];
        rdkit_vector_normalize(&mut data).unwrap();
        assert_eq!(data, [1.0]);
    }

    #[test]
    fn vector_normalize_source_negative_component() {
        let mut data = [-1.5e-16];
        rdkit_vector_normalize(&mut data).unwrap();
        assert_eq!(data, [-1.0]);
    }

    #[test]
    fn vector_normalize_source_nan_comparison_and_division() {
        let mut single = [f64::NAN];
        rdkit_vector_normalize(&mut single).unwrap();
        assert!(single[0].is_nan());
        let mut mixed = [f64::NAN, 1.0];
        rdkit_vector_normalize(&mut mixed).unwrap();
        assert!(mixed.iter().all(|value| value.is_nan()));
    }

    #[test]
    fn vector_normalize_source_infinity_division() {
        let mut data = [f64::INFINITY, 1.0];
        rdkit_vector_normalize(&mut data).unwrap();
        assert!(data[0].is_nan());
        assert_eq!(data[1].to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn vector_normalize_source_zero_preserves_signed_zero_on_error() {
        let mut data = [0.0, -0.0];
        let before = data.map(f64::to_bits);
        assert_eq!(
            rdkit_vector_normalize(&mut data),
            Err(ConformerError::CannotNormalizeZeroLengthVector)
        );
        assert_eq!(data.map(f64::to_bits), before);
    }

    #[test]
    fn vector_normalize_source_empty_is_structural_error() {
        assert_eq!(
            rdkit_vector_normalize(&mut []),
            Err(ConformerError::CannotNormalizeZeroLengthVector)
        );
    }

    #[test]
    fn vector_normalize_source_random_caller_propagates_error() {
        assert_eq!(
            rdkit_vector_set_to_random(0, 42),
            Err(ConformerError::CannotNormalizeZeroLengthVector)
        );
    }

    #[test]
    fn vector_normalize_source_generation_retains_numeric_error() {
        let error = crate::generation::GenerationError::from(
            ConformerError::CannotNormalizeZeroLengthVector,
        );
        assert!(matches!(
            error,
            crate::generation::GenerationError::Numeric(
                ConformerError::CannotNormalizeZeroLengthVector
            )
        ));
    }

    const RDKIT_RANDOM_MODULUS: u64 = 2_147_483_647;
    const RDKIT_RANDOM_MULTIPLIER: u64 = 48_271;
    fn assert_close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() < 1.0e-12,
            "actual={actual} expected={expected}"
        );
    }

    #[derive(Debug)]
    struct FixedDoubleRng {
        values: Vec<f64>,
        idx: usize,
    }

    impl FixedDoubleRng {
        fn new(values: Vec<f64>) -> Self {
            Self { values, idx: 0 }
        }
    }

    impl RdkitDoubleRng for FixedDoubleRng {
        fn next_unit_f64(&mut self) -> f64 {
            let value = self.values[self.idx];
            self.idx += 1;
            value
        }
    }

    fn rdkit_minstd_next_raw(state: &mut u64) -> u32 {
        *state = (*state * RDKIT_RANDOM_MULTIPLIER) % RDKIT_RANDOM_MODULUS;
        *state as u32
    }

    fn rdkit_minstd_next_unit(state: &mut u64) -> f64 {
        let raw = rdkit_minstd_next_raw(state) as f64;
        (raw - 1.0) / (RDKIT_RANDOM_MODULUS as f64 - 1.0)
    }

    fn rdkit_minstd_seed_state(seed: i32) -> u64 {
        let mut state = seed as u64 % RDKIT_RANDOM_MODULUS;
        if state == 0 {
            state = 1;
        }
        state
    }

    fn pick_random_dist_mat_bounds() -> BoundsMatrix {
        let mut mmat = BoundsMatrix::new(3).expect("matrix");
        mmat.set_lower(1, 0, 1.0).expect("set lower");
        mmat.set_upper(1, 0, 3.0).expect("set upper");
        mmat.set_lower(2, 0, 2.0).expect("set lower");
        mmat.set_upper(2, 0, 6.0).expect("set upper");
        mmat.set_lower(2, 1, 4.0).expect("set lower");
        mmat.set_upper(2, 1, 8.0).expect("set upper");
        mmat
    }

    #[test]
    fn pick_random_dist_mat_rng_overload_fills_lower_triangle_in_source_order() {
        let mmat = pick_random_dist_mat_bounds();
        let mut dist_mat = SymmMatrix::new(3);
        let mut rng = FixedDoubleRng::new(vec![0.25, 0.5, 0.75]);

        let largest = pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);

        assert_eq!(dist_mat.get_data(), &[0.0, 1.5, 0.0, 4.0, 7.0, 0.0]);
        assert_eq!(dist_mat.get_val(1, 0), 1.5);
        assert_eq!(dist_mat.get_val(0, 1), 1.5);
        assert_eq!(largest, 7.0);
    }

    #[test]
    fn pick_random_dist_mat_seed_overload_reseeds_rdkit_minstd_rng() {
        let mmat = pick_random_dist_mat_bounds();
        let mut first = SymmMatrix::new(3);
        let mut second = SymmMatrix::new(3);

        let first_largest = pick_random_dist_mat(&mmat, &mut first, 1);
        let second_largest = pick_random_dist_mat(&mmat, &mut second, 1);

        assert_eq!(first.get_data(), second.get_data());
        assert_eq!(first_largest, second_largest);

        let raw = 48_271.0;
        let rval = (raw - 1.0) / (2_147_483_647.0 - 1.0);
        assert_eq!(first.get_val(1, 0), 1.0 + rval * 2.0);
    }

    #[test]
    fn rdkit_minstd_rand_matches_boost_seed_modulus_boundary() {
        let mut direct = RdkitRandomEngine::from_seed(1);
        let mut wrapped = RdkitRandomEngine::from_seed(2_147_483_647);

        assert_eq!(direct.next_unit_f64(), wrapped.next_unit_f64());
    }

    #[test]
    fn embedder_conformer_seed_policy_matches_rdkit_non_sequential_and_sequential_modes() {
        assert_eq!(
            rdkit_embedder_conformer_seed(-1, 0, false).expect("original checked source seed"),
            -1
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(0, 0, false).expect("original checked source seed"),
            0
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(7, 0, false).expect("original checked source seed"),
            7
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(7, 1, false).expect("original checked source seed"),
            14
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(7, 2, false).expect("original checked source seed"),
            21
        );

        assert_eq!(
            rdkit_embedder_conformer_seed(0, 0, true).expect("original checked source seed"),
            1
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(7, 0, true).expect("original checked source seed"),
            8
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(7, 2, true).expect("original checked source seed"),
            10
        );
    }

    #[test]
    fn embedder_conformer_seed_policy_matches_rdkit_overflow_hash_branch() {
        assert!(rdkit_embedder_multiplication_overflows(46_342, 46_341));
        assert_eq!(
            rdkit_embedder_conformer_seed(46_341, 46_341, false)
                .expect("original checked source seed"),
            143_656
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(100_000, 100_000, false)
                .expect("original checked source seed"),
            1_410_365_413
        );
    }

    #[test]
    fn vector_set_to_random_clock_seeded_path_returns_normalized_vector() {
        let vec = rdkit_vector_set_to_random(3, 0).expect("clock-seeded vector");
        let norm_sq = vec.iter().map(|value| value * value).sum::<f64>();

        assert_eq!(vec.len(), 3);
        assert!(vec.iter().all(|value| value.is_finite()));
        assert!((norm_sq - 1.0).abs() < 1.0e-9);
    }

    #[test]
    fn power_eigen_solver_accepts_clock_seeded_path() {
        let mut mat = symm_matrix_from_distances(1, &[(0, 0, 1.0)]);
        let mut eigen_values = vec![0.0; 1];

        assert!(power_eigen_solver(1, &mut mat, &mut eigen_values, None, 0).expect("power solver"));
        assert!(eigen_values[0].is_finite());
    }

    #[test]
    fn pick_random_dist_mat_seed_overload_preserves_unseeded_global_stream() {
        let mmat = pick_random_dist_mat_bounds();
        let mut seeded = SymmMatrix::new(3);
        let mut continued = SymmMatrix::new(3);

        pick_random_dist_mat(&mmat, &mut seeded, 7);
        pick_random_dist_mat(&mmat, &mut continued, -1);

        assert_ne!(seeded.get_data(), continued.get_data());
    }

    #[test]
    #[should_panic(expected = "Size mismatch")]
    fn pick_random_dist_mat_rejects_size_mismatch() {
        let mmat = pick_random_dist_mat_bounds();
        let mut dist_mat = SymmMatrix::new(2);
        let mut rng = FixedDoubleRng::new(vec![0.0]);

        pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);
    }

    #[test]
    #[should_panic]
    fn pick_random_dist_mat_rejects_upper_bound_below_lower_bound() {
        let mut mmat = BoundsMatrix::new(2).expect("matrix");
        mmat.set_lower(1, 0, 3.0).expect("set lower");
        mmat.set_upper(1, 0, 2.0).expect("set upper");
        let mut dist_mat = SymmMatrix::new(2);
        let mut rng = FixedDoubleRng::new(vec![0.0]);

        pick_random_dist_mat_with_rng(&mmat, &mut dist_mat, &mut rng);
    }

    #[test]
    fn symm_matrix_set_val_uses_same_lower_triangular_storage_for_either_index_order() {
        let mut mat = SymmMatrix::new(3);

        mat.set_val(0, 2, 4.5);

        assert_eq!(mat.get_val(2, 0), 4.5);
        assert_eq!(mat.get_data()[3], 4.5);
    }

    fn symm_matrix_from_distances(n: usize, distances: &[(usize, usize, f64)]) -> SymmMatrix {
        let mut mat = SymmMatrix::new(n);
        for &(i, j, value) in distances {
            mat.set_val(i, j, value);
        }
        mat
    }

    fn point_distance(a: &[f64], b: &[f64]) -> f64 {
        a.iter()
            .zip(b)
            .map(|(ai, bi)| {
                let d = ai - bi;
                d * d
            })
            .sum::<f64>()
            .sqrt()
    }

    #[test]
    fn compute_initial_coords_embeds_two_point_distance() {
        let dist_mat = symm_matrix_from_distances(2, &[(1, 0, 2.0)]);
        let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];

        assert!(
            compute_initial_coords(&dist_mat, &mut positions, false, 2, 11)
                .expect("compute initial coords")
        );

        assert_close(point_distance(&positions[0], &positions[1]), 2.0);
    }

    #[test]
    fn compute_initial_coords_rng_overload_rejects_size_mismatch() {
        let dist_mat = symm_matrix_from_distances(2, &[(1, 0, 2.0)]);
        let mut positions = vec![vec![0.0; 3]];
        let mut rng = FixedDoubleRng::new(vec![0.25]);

        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            compute_initial_coords_with_rng(&dist_mat, &mut positions, &mut rng, false, 2)
                .expect("compute initial coords");
        }));

        assert!(result.is_err());
    }

    #[test]
    fn compute_initial_coords_fails_when_zero_eigen_threshold_is_reached_for_more_than_three_points()
     {
        let dist_mat = SymmMatrix::new(4);
        let mut positions = vec![vec![0.0; 3]; 4];

        assert!(
            !compute_initial_coords(&dist_mat, &mut positions, false, 2, 17)
                .expect("compute initial coords")
        );
    }

    #[test]
    fn compute_initial_coords_fails_negative_eigenvalue_without_randomization() {
        let dist_mat = symm_matrix_from_distances(3, &[(1, 0, 1.0), (2, 0, 1.0), (2, 1, 3.0)]);
        let mut positions = vec![vec![0.0; 3]; 3];

        assert!(
            !compute_initial_coords(&dist_mat, &mut positions, false, 2, 19)
                .expect("compute initial coords")
        );
    }

    #[test]
    fn compute_initial_coords_randomizes_negative_eigenvalue_when_requested() {
        let dist_mat = symm_matrix_from_distances(3, &[(1, 0, 1.0), (2, 0, 1.0), (2, 1, 3.0)]);
        let mut positions = vec![vec![0.0; 3]; 3];
        let mut rng = FixedDoubleRng::new(vec![0.25, 0.75, 0.5, 0.125, 0.875, 0.625]);

        assert!(
            compute_initial_coords_with_rng(&dist_mat, &mut positions, &mut rng, true, 2)
                .expect("compute initial coords")
        );

        assert!(positions.iter().flatten().all(|coord| coord.is_finite()));
    }

    #[test]
    fn compute_random_coords_rng_overload_fills_points_in_source_order_exactly() {
        let mut positions = vec![vec![0.0; 3], vec![0.0; 2], vec![0.0; 1]];
        let mut rng = FixedDoubleRng::new(vec![0.0, 0.25, 0.5, 0.75, 1.0, 0.125]);

        assert!(compute_random_coords_with_rng(
            &mut positions,
            4.0,
            &mut rng
        ));

        assert_eq!(
            positions,
            vec![vec![-2.0, -1.0, 0.0], vec![1.0, 2.0], vec![-1.5]]
        );
        assert_eq!(rng.idx, 6);
    }

    #[test]
    fn compute_random_coords_seed_overload_reseeds_boost_minstd_rng_exactly() {
        let mut positions = vec![vec![0.0; 3], vec![0.0; 3]];
        let mut expected_state = rdkit_minstd_seed_state(42);
        let expected = (0..2)
            .map(|_| {
                (0..3)
                    .map(|_| 6.0 * (rdkit_minstd_next_unit(&mut expected_state) - 0.5))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        assert!(compute_random_coords(&mut positions, 6.0, 42));

        assert_eq!(positions, expected);
    }

    #[test]
    fn compute_random_coords_seed_boundary_matches_boost_seed_zero_adjustment() {
        let mut wrapped_positions = vec![vec![0.0; 2]];
        let mut direct_positions = vec![vec![0.0; 2]];

        assert!(compute_random_coords(
            &mut wrapped_positions,
            2.0,
            2_147_483_647
        ));
        assert!(compute_random_coords(&mut direct_positions, 2.0, 1));

        assert_eq!(wrapped_positions, direct_positions);
    }

    #[test]
    fn compute_random_coords_empty_positions_returns_true_without_consuming_rng() {
        let mut positions: Vec<Vec<f64>> = Vec::new();
        let mut rng = FixedDoubleRng::new(vec![0.25]);

        assert!(compute_random_coords_with_rng(
            &mut positions,
            1.0,
            &mut rng
        ));

        assert!(positions.is_empty());
        assert_eq!(rng.idx, 0);
    }

    #[test]
    #[should_panic(expected = "bad boxSize")]
    fn compute_random_coords_rejects_non_positive_box_size() {
        let mut positions = vec![vec![0.0; 3]];
        let mut rng = FixedDoubleRng::new(vec![0.25]);

        compute_random_coords_with_rng(&mut positions, 0.0, &mut rng);
    }
    #[test]
    fn compute_initial_coords_power_solver_follows_source_seeded_diagonal_iteration() {
        let mut mat = symm_matrix_from_distances(3, &[(0, 0, 5.0), (1, 1, 3.0), (2, 2, 1.0)]);
        let mut eigen_values = vec![0.0; 2];
        let mut eigen_vectors = DoubleMatrix::new(2, 3);

        assert!(
            power_eigen_solver(2, &mut mat, &mut eigen_values, Some(&mut eigen_vectors), 23)
                .expect("power solver")
        );

        assert!((eigen_values[0] - 3.0).abs() < 1.0e-3);
        assert!((eigen_values[1] - 1.006_175_852_283_403_2).abs() < 1.0e-3);
    }

    #[test]
    fn embedder_conformer_seed_source_release_unsigned_wrap_and_invariants() {
        assert_eq!(
            rdkit_embedder_conformer_seed(-1, usize::MAX, false).unwrap(),
            -1
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(0, usize::MAX, false).unwrap(),
            0
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(0, u32::MAX as usize - 1, true).unwrap(),
            -1
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(0, u32::MAX as usize, true).unwrap(),
            0
        );
        let error = rdkit_embedder_conformer_seed(i32::MAX, 0, true).unwrap_err();
        assert_eq!(
            error,
            EmbedSeedError::ResultInvariant {
                seed: i32::MAX,
                index: 0,
                value: i32::MIN
            }
        );
        assert_eq!(
            error.to_string(),
            "Something went wrong calculating a new seed"
        );
        assert_eq!(
            rdkit_embedder_conformer_seed(-2, 0, false).unwrap_err(),
            EmbedSeedError::InvalidRandomSeed { seed: -2 }
        );
    }
}
