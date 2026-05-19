# DG Bounds Source-Closure Re-Audit

Date: 2026-05-18

Baseline reviewed against:
- `third_party/rdkit/Code/DistGeom/BoundsMatrix.h`
- `third_party/rdkit/Code/DistGeom/TriangleSmooth.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp::getMolBoundsMatrix(...)`

## Scope verdict

The redesigned Rust DG bounds implementation has substantial coverage, but the
selected source-closure slice through planned Step 50 is still open. The
remaining direct gaps inside that slice are concrete and local, so the existing
plan remains executable without regeneration.

## Function-by-function gaps relevant to Steps 5-50

### BoundsMatrix.h

- `BoundsMatrix::setUpperBoundIfBetter(...)`
  - Missing as a source-anchored helper on `BoundsMatrix`.
  - Current callers only have `set_upper(...)`; no exact RDKit conditional
    setter exists.
- `BoundsMatrix::setLowerBoundIfBetter(...)`
  - Missing as a source-anchored helper on `BoundsMatrix`.
- `BoundsMatrix::checkValid()`
  - Missing as a source-anchored helper on `BoundsMatrix`.

### TriangleSmooth.cpp

- `triangleSmoothBounds(BoundsMatPtr, double)`
  - Missing forwarding overload.
  - Current Rust code has only the pointer-like in-place smoother method on
    `BoundsMatrix`.
- `triangleSmoothBounds(BoundsMatrix *, double)`
  - Present as `BoundsMatrix::triangle_smooth(...)`.
  - Behavior looks aligned, but Step 11 requires pointer/pointer-like entrypoint
    alignment with the RDKit overload pair.

### AtomTyper.cpp helpers used by DG bounds

- `UFF::Tools::addAtomChargeFlags(...)`
  - Current Rust function is summary-style and not source-framed line-by-line.
  - It omits copied branch comments and exact charge-mismatch/toleration shape.
  - Lanthanide and special-element labeling needs direct-source re-expression.
- `UFF::Tools::getAtomLabel(...)`
  - Current Rust function is summary-style and not source-framed line-by-line.
  - Hybridization/aromatic/conjugation label suffix logic needs exact source
    framing.

### Geometry helper

- `RDGeom::compute13Dist`
  - Current Rust helper is only a summary comment plus formula.
  - It needs the exact source frame copied inline and targeted tests at this
    helper boundary.

### Topological distance helper used by 1-4 / 1-5

- Current Rust all-pairs helper uses Floyd-Warshall over an adjacency matrix.
- RDKit path-distance input is conceptually `MolOps::getDistanceMat(mol)`.
- For this plan slice, the gap is not missing functionality but missing direct
  source anchoring and helper closure for the shortest-topological-distance
  behavior that feeds `set14Bounds(...)` and `set15Bounds(...)`.

### `set12Bounds(...)`

- Current function already contains a more complete port than the summary block
  above it, but still carries `RDKit✔️❗` markers.
- Remaining closure work:
  - tighten source reproduction frames around UFF typing / fallback branches;
  - verify the conjugated-5-ring squish branch and visited12 bookkeeping with
    targeted tests;
  - upgrade or preserve markers with evidence instead of mixed summary markers.

### `set13Bounds(...)`

- Current function still carries a top-level `RDKit✔️❌` block.
- The body is implemented, but the remaining gap is source-closure:
  - copied source frame is still summarized with `{ ... }` elisions;
  - marker state overstates unclosed source reproduction;
  - targeted tests for newly closed branches are still missing.

## Out-of-slice observations

- Large 1-4 / 1-5 helpers after Step 50 still retain many `RDKit✔️❌` /
  `RDKit❗✔️` markers and are intentionally deferred to later plan steps.
- `collectBondsAndAngles(...)` and `setTopolBounds(...)` are already closer to
  closure than the Step 5-50 targets.

## Recommended execution order

Execute the existing plan as written through Step 50. No regeneration is needed
for the first 50 steps because each remaining gap is bounded to a single helper
or a single already-implemented function whose source framing and tests can be
closed incrementally.
