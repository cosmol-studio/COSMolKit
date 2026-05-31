# RDKit Conformer Generation Final Marker Audit

Date: 2026-05-30

Scope audited:

- `crates/cosmolkit-core/src/chemistry/distgeom.rs`
- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_preferences.rs`
- `crates/cosmolkit-core/tests/rdkit_conformer_generation_parity.rs`

Marker counts from the audited conformer-generation surface:

- `RDKit❌❌`: 0
- `RDKit❗❗`: 0
- `RDKit❗✔️`: 65
- `RDKit✔️❌`: 490
- `RDKit✔️❗`: 98

## 1. First-axis behavioral gap closure

No first-axis `RDKit❌❌` block remains in the audited conformer-generation
surface.

The previously remaining first-axis gap was RDKit's
`MolAlign::details::symmetrizeTerminalAtoms()` branch used by symmetry-aware
pruning when `symmetrizeConjugatedTerminalGroupsForPruning == true`.

Current COSMolKit status:

- The branch is now source-ported in
  `crates/cosmolkit-core/src/chemistry/distgeom.rs`.
- The helper mutates the query-side clone exactly along the RDKit source shape:
  neutralize terminal formal charge and replace the matched terminal bond with a
  single-or-double query bond before self-matching.
- Local regression coverage now asserts the positive behavior at
  `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs`.
- Multi-conformer parity coverage no longer expects this branch to fail closed.

Assessment:

- No first-axis `RDKit❌❌` block remains in the audited conformer-generation
  path.

## 2. First-axis gap closed during this pass

The previous unseeded `clock()`-based RNG paths are no longer `RDKit❌❌`.

Closed paths:

- `RDNumeric::Vector::setToRandom` clock seeding in
  `crates/cosmolkit-core/src/chemistry/distgeom.rs:5823`
- `RDNumeric::EigenSolvers::powerEigenSolver` clock seeding in
  `crates/cosmolkit-core/src/chemistry/distgeom.rs:5906`

Current status:

- Both paths now call the C runtime `clock()` through a Rust `extern "C"`
  binding and remain source-backed.
- Regression tests were updated to assert acceptance instead of panic in
  `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs:3033` and
  `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs:3043`.

## 3. Remaining `RDKit❗✔️` behavioral-partial markers

No `RDKit❗❗` block remains, but `RDKit❗✔️` markers still indicate behavior that
is not yet claimed as fully reproduced even though local inspection supports the
current complexity shape.

The remaining `RDKit❗✔️` markers are concentrated in these source regions:

- `crates/cosmolkit-core/src/chemistry/distgeom.rs:43`
  `BoundsMatrixBuilder.cpp` source import boundary.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:4569`
  VDW radii provenance note.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:5048`
  ideal bond-angle derivation and `_setRingAngle`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:7593`
  `set13BoundsHelper`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:7679`
  `Path14Configuration`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:7689`
  `Path14Type`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:7704`
  `ComputedData`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:9652`
  `set15Bounds`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:9835`
  detailed `set15Bounds` branch notes.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:9979`
  `set14Bounds`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:10268`
  `setLowerBoundVDW`.

Assessment:

- These are the remaining behavior-not-yet-fully-upgraded marker groups in the
  bounds-builder subgraph.
- They are not fail-open unsupported paths, but they are also not yet
  marker-closed first-axis parity claims.

## 4. Remaining `RDKit✔️❗` performance-review-unresolved markers

The unresolved second-axis markers are concentrated in two places:

- `crates/cosmolkit-core/src/chemistry/distgeom.rs`
  large bounds-builder regions including `set12Bounds` and `set13Bounds`.
- `crates/cosmolkit-core/src/chemistry/forcefield/crystalff/torsion_preferences.rs:308`
  `getExperimentalTorsions`.

Assessment:

- These markers do not assert a behavioral gap.
- They do assert that local complexity/performance review is not closed.
- They should not be collapsed into blanket parity language in public support
  docs.

## 5. Remaining `RDKit✔️❌` known complexity/performance gaps

The known second-axis regressions are concentrated in the bounds-matrix helper
surface:

- `crates/cosmolkit-core/src/chemistry/distgeom.rs:7191`
  `triangleSmoothBounds`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:8155`
  `_record14Path`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:9241`
  `_checkH2NX3H1OX2`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:9272`
  `_checkNhChChNh`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:9328`
  `_checkMacrocycleAllInSameRingAmideEster14`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:9458`
  `_checkAmideEster15`.
- `crates/cosmolkit-core/src/chemistry/distgeom.rs:10495`
  `getMolBoundsMatrix(...)` wrapper/export path.

Assessment:

- These are second-axis-only gaps.
- They do not currently invalidate deterministic seeded single-conformer parity
  or deterministic batch seed policy.

## 6. Summary

Current audited conformer-generation state after this pass:

- Seeded RNG path: source-backed, deterministic, and no remaining `clock()`
  first-axis gap.
- Unseeded `clock()` path: source-backed, no longer panic-based.
- Symmetry pruning with terminal-group symmetrization: source-ported and no
  longer a first-axis gap.
- Bounds-builder helper surface: still contains unresolved `RDKit❗✔️`,
  `RDKit✔️❗`, and `RDKit✔️❌` markers, so the feature cannot be described as
  fully marker-closed.
