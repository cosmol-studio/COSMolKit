# DG Bounds RDKit Remaining Source Scan

Date: 2026-05-17

## Purpose

This report re-audits the DG bounds port after checklist steps `1-139`.

The audit baseline is the DG bounds call graph intentionally tracked by
`dev/dg_bounds_rdkit_full_port_checklist.md`:

- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`
- `third_party/rdkit/Code/DistGeom/BoundsMatrix.h`
- `third_party/rdkit/Code/DistGeom/TriangleSmooth.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp`

The audited Rust implementation lives in:

- `crates/cosmolkit-core/src/chemistry/distgeom.rs`

## Current Verdict

Within the selected DG bounds scope, remaining direct source gaps are now
**zero**.

That verdict is limited to the baseline defined by the checklist:

1. DG bounds matrix initialization and bound-setting helpers
2. DG 1-2 / 1-3 / 1-4 / 1-5 topology-bound generation
3. raw `BoundsMatrix` upper/lower triangle semantics required by that call graph
4. triangle smoothing semantics required by that call graph
5. the public `getMolBoundsMatrix(...)` wrapper behavior mapped to
   `dg_bounds_matrix_with_options()` / `dg_bounds_matrix()`

This report does **not** make any RDKit parity claim beyond source-backed port
closure plus local targeted tests.

## Source-Closed Coverage

The following direct RDKit functions are present in Rust with copied source
anchors and corresponding behavior/tests in the current DG bounds scope:

### `BoundsMatrix.h`

1. `DistGeom::BoundsMatrix::getUpperBound(...)`
2. `DistGeom::BoundsMatrix::setUpperBound(...)`
3. `DistGeom::BoundsMatrix::setLowerBound(...)`
4. `DistGeom::BoundsMatrix::getLowerBound(...)`

### `TriangleSmooth.cpp`

5. `DistGeom::triangleSmoothBounds(...)` core `BoundsMatrix *` behavior

### `BoundsMatrixBuilder.cpp`

6. `DGeomHelpers::ComputedData::visitedBound(...)`
7. `DGeomHelpers::initBoundsMat(DistGeom::BoundsMatrix *mmat, ...)`
8. `DGeomHelpers::initBoundsMat(DistGeom::BoundsMatPtr mmat, ...)`
9. `DGeomHelpers::_checkAndSetBounds(...)`
10. `DGeomHelpers::set12Bounds(...)`
11. `DGeomHelpers::isHBondAcceptor(...)`
12. `DGeomHelpers::isHBondDonor(...)`
13. `DGeomHelpers::isHinHBondDonor(...)`
14. `DGeomHelpers::setLowerBoundVDW(...)`
15. `DGeomHelpers::isLargerSP2Atom(...)`
16. `DGeomHelpers::_set13BoundsHelper(...)`
17. `DGeomHelpers::_setRingAngle(...)`
18. `DGeomHelpers::set13Bounds(...)`
19. `DGeomHelpers::_record14Path(...)`
20. `DGeomHelpers::_setInRing14Bounds(...)`
21. `DGeomHelpers::_setTwoInSameRing14Bounds(...)`
22. `DGeomHelpers::_setTwoInDiffRing14Bounds(...)`
23. `DGeomHelpers::_setShareRingBond14Bounds(...)`
24. `DGeomHelpers::_checkH2NX3H1OX2(...)`
25. `DGeomHelpers::_checkNhChChNh(...)`
26. `DGeomHelpers::_checkAmideEster14(...)`
27. `DGeomHelpers::_checkMacrocycleAllInSameRingAmideEster14(...)`
28. `DGeomHelpers::_isCarbonyl(...)`
29. `DGeomHelpers::_checkAmideEster15(...)`
30. `DGeomHelpers::_setChain14Bounds(...)`
31. `DGeomHelpers::_checkMacrocycleTwoInSameRingAmideEster14(...)`
32. `DGeomHelpers::_setMacrocycleTwoInSameRing14Bounds(...)`
33. `DGeomHelpers::_setMacrocycleAllInSameRing14Bounds(...)`
34. `DGeomHelpers::set14Bounds(...)`
35. `DGeomHelpers::_compute15DistsCisCis(...)`
36. `DGeomHelpers::_compute15DistsCisTrans(...)`
37. `DGeomHelpers::_compute15DistsTransTrans(...)`
38. `DGeomHelpers::_compute15DistsTransCis(...)`
39. `DGeomHelpers::_set15BoundsHelper(...)`
40. `DGeomHelpers::set15Bounds(...)`
41. `DGeomHelpers::collectBondsAndAngles(...)`
42. `DGeomHelpers::setTopolBounds(...)` first overload
43. `DGeomHelpers::setTopolBounds(...)` second overload

### `rdDistGeom.cpp`

44. `RDKit::getMolBoundsMatrix(...)`

## Resolved Historical Blockers

The previous blocking gaps are now closed:

1. `BoundsMatrix` no longer uses split symmetric `lower/upper` storage.
   The Rust implementation now stores the same raw upper-triangle/lower-triangle
   layout required by RDKit semantics.
2. `triangleSmoothBounds(...)` no longer uses the earlier simplified
   multi-pass approximation. The current Rust implementation follows the RDKit
   `k/i/j` traversal, the two lower-difference formulas, the tolerance
   reconciliation branch, and the failure return.
3. `dg_bounds_matrix()` no longer uses the old wrapper defaults. The current
   default call matches the RDKit wrapper settings in scope:
   `set15bounds=true`, `scale_vdw=false`, `do_triangle_smoothing=true`,
   `use_macrocycle14config=false`.
4. Legacy local `distgeom` tests that still assumed the pre-port matrix/export
   semantics were updated to assert the current source-backed DG bounds path.

## Scope Boundary Notes

The following RDKit functions are **not** treated as remaining DG bounds gaps in
this report because they are outside the selected baseline or are convenience
wrappers that do not add new DG bounds behavior beyond the already ported core:

1. `DistGeom::BoundsMatrix::setUpperBoundIfBetter(...)`
2. `DistGeom::BoundsMatrix::setLowerBoundIfBetter(...)`
3. `DistGeom::BoundsMatrix::checkValid()`
4. `DistGeom::triangleSmoothBounds(BoundsMatPtr, ...)` forwarding overload
5. Python-binding glue in `Wrap/rdDistGeom.cpp` unrelated to DG bounds
   generation semantics

If the project later expands the DG bounds baseline to include those APIs as
first-class Rust surfaces, this report must be regenerated.

## Validation Used For This Rescan

The local targeted validation relevant to this audit is:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict distgeom -- --nocapture
```

Result at audit time:

- `66 passed; 0 failed`

## Zero-Gap Status

For the DG bounds baseline tracked by the checklist, zero-gap status is
**granted**.

That means:

1. every baseline direct helper/function has a source-backed Rust
   implementation,
2. the matrix-layer and smoothing-layer semantics that previously blocked
   closure are now ported,
3. the local strict `distgeom` test slice passes after the legacy test repair.
