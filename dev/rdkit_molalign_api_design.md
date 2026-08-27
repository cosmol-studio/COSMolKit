# MolAlign / RMSD API Design

## Boundary

This design covers ordinary RDKit `MolAlign` and `AlignPoints` behavior only.
O3A, `RandomTransform`, and MMFF/Crippen-based Open3DAlign are separate future
capabilities. The implementation uses one private source-backed Rust numerical
core and does not add an FFI dependency because the pinned `AlignPoints.cpp`
algorithm is directly reproducible with Rust scalar operations.

## Core Types

```rust
pub struct AlignmentAtomMap {
    pub probe_atom: usize,
    pub reference_atom: usize,
}

pub struct AlignmentParameters {
    pub probe_conformer_id: i32,      // -1 selects first stored conformer
    pub reference_conformer_id: i32,  // -1 selects first stored conformer
    pub atom_map: Option<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
}

pub struct BestAlignmentParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub reflect: bool,
    pub max_iterations: u32,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
    pub ignore_hydrogens: bool,
    pub num_threads: i32,
}

pub struct AllConformerRmsdParameters {
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
    pub ignore_hydrogens: bool,
    pub num_threads: i32,
}

pub struct CoordinateRmsdParameters {
    pub probe_conformer_id: i32,
    pub reference_conformer_id: i32,
    pub atom_maps: Vec<Vec<AlignmentAtomMap>>,
    pub weights: Option<Vec<f64>>,
    pub max_matches: i32,
    pub symmetrize_conjugated_terminal_groups: bool,
}

pub struct AlignmentTransform {
    pub matrix: [[f64; 4]; 4],
}

pub struct AlignmentResult {
    pub rmsd: f64,
    pub transform: AlignmentTransform,
    pub atom_map: Vec<AlignmentAtomMap>,
}

pub struct ConformerRmsd {
    pub probe_conformer_id: usize,
    pub reference_conformer_id: usize,
    pub rmsd: f64,
}

pub struct ConformerAlignmentReport {
    pub rmsds: Vec<f64>,
}
```

The exact field names may be adjusted to the existing public naming style, but
the semantic distinctions are fixed: an alignment result contains the chosen
map and transform; a coordinate-frame RMS result contains only a measurement;
and conformer-pair output retains real stored IDs and source ordering.

## Read-Only Methods

The following methods never mutate either receiver, including coordinates,
properties, topology, caches, or conformer metadata:

```rust
Molecule::alignment_transform_to(&self, reference: &Molecule, params: &AlignmentParameters)
Molecule::coordinate_rmsd_to(&self, reference: &Molecule, params: &CoordinateRmsdParameters)
Molecule::best_alignment_to(&self, reference: &Molecule, params: &BestAlignmentParameters)
Molecule::best_rmsd_to(&self, reference: &Molecule, params: &BestAlignmentParameters)
Molecule::all_conformer_best_rmsds(&self, params: &AllConformerRmsdParameters)
```

`coordinate_rmsd_to` does not apply a transform. `alignment_transform_to` and
`best_alignment_to` compute transforms but do not apply them. Automatic mapping
uses the existing matcher with source-specific parameters. The all-conformer
parameter type deliberately omits conformer IDs, reflection, and iteration
limits because pinned `getAllConformerBestRMS` fixes those choices internally.

## Explicit Coordinate Mutation

Mutation is available only through clearly named methods and registered
coordinate operations:

```rust
Molecule::with_alignment_to(&self, reference: &Molecule, params: &AlignmentParameters)
Molecule::align_to_(&mut self, reference: &Molecule, params: &AlignmentParameters)
Molecule::with_aligned_conformers_with_params(&self, params: ConformerAlignmentParameters)
    -> Result<(Molecule, ConformerAlignmentReport), OperationError>
Molecule::align_conformers_with_params_(&mut self, params: ConformerAlignmentParameters)
    -> Result<ConformerAlignmentReport, OperationError>
```

Value-style methods return a new molecule plus the alignment result. In-place
methods end in `_` and return the measurement/result after successful
coordinate mutation. The report is the ordered RDKit `RMSlist` produced by the
same registered coordinate transaction; it is neither recalculated nor exposed
if operation-contract or molecule-invariant validation fails. No method named
`rmsd`, `calculate`, or `measure` mutates.

## Validation And Errors

The structured error covers missing/invalid conformer IDs, empty point maps,
atom-map bounds, map/weight length mismatch, non-positive weights, missing
automatic matches, worker termination, and coordinate-row length mismatches.
`num_threads` follows RDKit's signed `getNumThreadsToUse()` semantics, including
zero and negative hardware-relative values. `max_matches` retains the source
`int` storage and legacy unsigned-conversion behavior. RDKit source
preconditions are represented as errors at the public boundary; they are not
panics, silent empty results, or fallback maps.

## Lifecycle

Alignment reads topology and coordinates only. Coordinate mutation invalidates
drawing state according to the existing coordinate-operation policy and
preserves topology, properties, and CIP state. Read-only operations do not
alter caches. A symmetry clone is temporary and never becomes user-visible
molecule state.

## Numerical Boundary

There is one implementation: the source-ordered pure Rust implementation.
Pinned RDKit expected data is generated outside the production dependency graph
and compares matrix, SSR/RMSD, map, ordering, and errors. FFI may only be added
for a specific numerical path after parity measurement proves that Rust cannot
reproduce the required result; source availability, test convenience, and a
generic dual-backend design are not sufficient reasons. That condition is not
present in the audited ordinary MolAlign boundary.

## Validation Status

The design is implemented through one pure-Rust numerical core and one public
Rust/Python API path. Focused, 152-row, 5,000-row, and complete ChEMBL 37 gates
pass. The complete audit processes 2,854,362 eligible molecules and compares
11,417,207 RMSD, transform, map, coordinate, order, and immutability fields
with zero mismatch. The source and architecture closure is recorded in
`dev/gap_reports/rdkit_molalign_full_port_validation.md`.
