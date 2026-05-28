# RDKit ForceField Core Symbol Inventory

Source files audited:

- `third_party/rdkit/Code/ForceField/ForceField.h`
- `third_party/rdkit/Code/ForceField/ForceField.cpp`
- `third_party/rdkit/Code/ForceField/Contrib.h`

## Namespaces And Free Functions

- `RDKit::ForceFieldsHelper::normalizeAngleDeg(double&)`
  - Branches: `angleDeg < -180.0`, `angleDeg > 180.0`, in-range after `fmod`.
- `RDKit::ForceFieldsHelper::computeDihedral(const RDGeom::PointPtrVect&, unsigned int, unsigned int, unsigned int, unsigned int, double*, double*, RDGeom::Point3D[4], RDGeom::Point3D[2], double[2])`
  - Branches: delegates through four point pointers and optional output pointers.
- `RDKit::ForceFieldsHelper::computeDihedral(const double*, unsigned int, unsigned int, unsigned int, unsigned int, double*, double*, RDGeom::Point3D[4], RDGeom::Point3D[2], double[2])`
  - Branches: builds four `RDGeom::Point3D` values from flat coordinate storage and delegates.
- `RDKit::ForceFieldsHelper::computeDihedral(const RDGeom::Point3D*, const RDGeom::Point3D*, const RDGeom::Point3D*, const RDGeom::Point3D*, double*, double*, RDGeom::Point3D[4], RDGeom::Point3D[2], double[2])`
  - Preconditions: all four point pointers must be non-null.
  - Branches: optional local `r`, `t`, `d`, and `cosPhi` outputs; clamped cross-product lengths; clamped cosine; optional signed dihedral output.
- `ForceFieldsHelper::calcEnergy`
  - Callable wrapper around `ForceFields::ForceField::calcEnergy(double*)`.
- `ForceFieldsHelper::calcGradient`
  - Branches: zero gradient first; call `calcGrad`; scale all gradients by `0.1`; compute max gradient; optional repeated `0.5` scaling while `maxGrad * gradScale > 10.0`.

## Classes And Fields

- `ForceFields::ForceFieldContrib`
  - Constructors: default and owner-bound constructor.
  - Virtual destructor.
  - Pure virtual methods: `getEnergy(double*) const`, `getGrad(double*, double*) const`, `copy() const`.
  - Field: `ForceField *dp_forceField`.
  - Friend: `ForceField`.
- `ForceFields::ForceField`
  - Constructor: `ForceField(unsigned int dimension = 3)`.
  - Destructor: clears point/contrib storage and deletes distance matrix.
  - Copy constructor: copies dimension and point count, resets initialization and distance matrix, deep-copies contribs, rebinds copied contribs to the new owner.
  - Fields:
    - `unsigned int d_dimension`
    - `bool df_init`
    - `unsigned int d_numPoints`
    - `double *dp_distMat`
    - `RDGeom::PointPtrVect d_positions`
    - `ContribPtrVect d_contribs`
    - `INT_VECT d_fixedPoints`
    - `unsigned int d_matSize`

## ForceField Methods And Branches

- `initialize()`
  - Clears old distance matrix, stores point count, allocates triangular matrix, initializes matrix to `-1.0`, then marks initialized.
- `minimize(unsigned int maxIts, double forceTol, double energyTol)`
  - Delegates to snapshot overload with `snapshotFreq = 0` and null snapshots.
- `minimize(unsigned int snapshotFreq, RDKit::SnapshotVect*, unsigned int maxIts, double forceTol, double energyTol)`
  - Preconditions: initialized and point-count/position-size match.
  - Branches: empty contrib list returns success; otherwise scatters coordinates, runs `BFGSOpt::minimize`, gathers optimized coordinates, returns optimizer status.
- `calcEnergy(std::vector<double>*) const`
  - Preconditions: initialized.
  - Branches: empty contrib list returns zero; optional contrib vector is cleared/reserved and populated; positions are scattered into flat storage before contrib iteration.
- `calcEnergy(double*)`
  - Preconditions: initialized and non-null position vector.
  - Branches: distance matrix reset first; empty contrib list returns zero; otherwise sums contrib energies over provided coordinates.
- `calcGrad(double*) const`
  - Preconditions: initialized and non-null gradient vector.
  - Branches: empty contrib list returns without modifying gradient; otherwise scatter positions, call contrib gradients, zero fixed-point gradient ranges.
- `calcGrad(double*, double*)`
  - Preconditions: initialized and non-null position/gradient vectors.
  - Branches: empty contrib list returns without modifying gradient; otherwise call contrib gradients over provided coordinates and zero fixed-point gradient ranges.
- `distance(unsigned int, unsigned int, double*)`
  - Preconditions: initialized and both indices in range.
  - Branches: swaps indices for triangular storage when `j < i`; invariant-checks matrix index; lazy-computes only when cached value is negative; computes from stored positions or provided flat coordinates.
- `distance(unsigned int, unsigned int, double*) const`
  - Branches: delegates to `sqrt(distance2(...))` without updating cache.
- `distance2(unsigned int, unsigned int, double*) const`
  - Preconditions: initialized and both indices in range.
  - Branches: swaps indices when `j < i`; computes from stored positions or provided flat coordinates without updating cache.
- `positions()`
  - Mutable and const reference accessors.
- `contribs()`
  - Mutable and const reference accessors.
- `dimension()`
  - Returns `d_dimension`.
- `numPoints()`
  - Returns `d_numPoints`.
- `fixedPoints()`
  - Mutable and const reference accessors.
- `scatter(double*) const`
  - Preconditions: initialized and non-null position vector.
  - Branches: copies every coordinate by dimension and postcondition-checks final offset.
- `gather(double*)`
  - Preconditions: initialized and non-null position vector.
  - Branches: writes every coordinate by dimension back into stored points.
- `initDistanceMatrix()`
  - Preconditions: nonzero point count, allocated distance matrix, matrix capacity sufficient.
  - Branches: fills triangular matrix with `-1.0`.

## Source Test Fixtures

- `third_party/rdkit/Code/ForceField/catch_tests.cpp`
  - Exercises force-field dihedral helper behavior and contrib-facing geometric calculations.
- `third_party/rdkit/Code/Numerics/Optimizer/testOptimizer.cpp`
  - Exercises the BFGS optimizer used by `ForceField::minimize`.
- Constraint-specific and UFF/MMFF fixture files are outside the force-field core steps and belong to later plan sections.
