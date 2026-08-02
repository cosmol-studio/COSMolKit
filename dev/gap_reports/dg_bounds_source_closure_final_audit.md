# DG Bounds Source Closure Final Audit

Date: 2026-05-18

## Scope

This audit covers the DG bounds generation baseline used by
`dev/archive/plans/dg_bounds_rdkit_full_port_checklist.md` after Step 159 refreshed the
direct baseline inventory to include `DGeomHelpers::_getAtomStereo(...)`:

- `third_party/rdkit/Code/DistGeom/BoundsMatrix.h`
- `third_party/rdkit/Code/DistGeom/TriangleSmooth.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp::getMolBoundsMatrix(...)`

The audited Rust implementation is:

- `crates/cosmolkit-core/src/chemistry/distgeom.rs`

## Final Verdict

Within the selected DG bounds baseline, this port no longer has any remaining
first-axis `RDKit❌*` markers in the active call chain. The remaining markers
are all one of:

1. `RDKit✔️❌`
   Behavior is treated as reproduced for the currently modeled input state
   space, but local review found a clear performance or data-structure cost.
2. `RDKit✔️❗`
   Behavior is treated as reproduced, but local performance equivalence was not
   granted after inspection.
3. `RDKit❗✔️`
   The remaining uncertainty is confined to local helper abstraction or source
   framing, not to an identified missing DG bounds branch in the selected
   baseline.

So the exact final state is:

- no known silent behavioral hole in the selected DG bounds call chain,
- no claim of zero marker residue,
- no claim of full RDKit parity outside the audited baseline,
- several deliberate performance/abstraction markers still remain and must stay
  visible.

## Function-by-Function Verdict

### Direct RDKit baseline functions still carrying `RDKit✔️❌`

1. `DistGeom::triangleSmoothBounds(BoundsMatrix *, double)`
   - Rust surface: `BoundsMatrix::triangle_smooth(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: behavior is source-backed and tested, but the Rust `Vec<Vec<f64>>`
     storage and checked indexing are an obvious cost relative to RDKit's flat
     contiguous matrix buffer.

2. `DGeomHelpers::_record14Path(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: behavior is covered, but path-flag recording uses `Vec<u64>` plus
     linear `contains()` checks where RDKit uses set/bitset-like storage.

3. `DGeomHelpers::_setInRing14Bounds(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: behavior is source-backed, but repeated helper lookups, `Vec`
     path-flag scans, and Rust-side branching/allocation shape are materially
     worse than the original.

4. `DGeomHelpers::_setTwoInSameRing14Bounds(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: same performance cost class as `_setInRing14Bounds(...)`.

5. `DGeomHelpers::_setMacrocycleTwoInSameRing14Bounds(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: same repeated lookup / Vec-backed path bookkeeping cost.

6. `DGeomHelpers::_setMacrocycleAllInSameRing14Bounds(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: same repeated lookup / Vec-backed path bookkeeping cost.

7. `DGeomHelpers::_setChain14Bounds(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: same repeated lookup / Vec-backed path bookkeeping cost.

8. `DGeomHelpers::_checkH2NX3H1OX2(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: behavior is tiny and covered, but the Rust atom-state access path
     is not a clear complexity match to RDKit's direct atom methods.

9. `DGeomHelpers::_checkNhChChNh(...)`
   - Marker: `RDKit✔️❌`
   - Verdict: keep as-is.
   - Reason: helper behavior is source-backed; the second axis stays negative
     because it is composed through the slower local helper path above.

10. `DGeomHelpers::_checkMacrocycleAllInSameRingAmideEster14(...)`
    - Marker: `RDKit✔️❌`
    - Verdict: keep as-is.
    - Reason: neighbor scans and bond lookups are materially less direct than
      RDKit's iterator-level access.

11. `DGeomHelpers::_checkAmideEster15(...)`
    - Marker: `RDKit✔️❌`
    - Verdict: keep as-is.
    - Reason: same helper-scan cost; behavior is covered.

12. `RDKit::getMolBoundsMatrix(...)`
    - Rust surface: `dg_bounds_matrix_with_options(...)`
    - Marker: `RDKit✔️❌`
    - Verdict: keep as-is.
    - Reason: the wrapper preserves control flow and defaults, but exporting as
      nested `Vec<Vec<f64>>` is materially less efficient than RDKit's direct
      contiguous NumPy allocation plus `memcpy()`.

### Direct RDKit baseline functions still carrying `RDKit✔️❗`

1. `DGeomHelpers::set12Bounds(...)`
   - Marker: `RDKit✔️❗`
   - Verdict: keep as-is.
   - Reason: behavior is reproduced for the modeled state space, but local
     review did not justify claiming performance equivalence across the UFF
     typing path, fallback branch, and squish bookkeeping.

2. `DGeomHelpers::set13Bounds(...)`
   - Marker: `RDKit✔️❗`
   - Verdict: keep as-is.
   - Reason: behavior is reproduced and tested, but the Rust port still does
     repeated `bond_between_idx_simple()` lookups inside nested loops where
     RDKit walks bond iterators and matrix cells more directly.

### Remaining `RDKit❗✔️` helpers in the DG bounds call chain

These are not treated as known missing baseline behaviors, but the first-axis
marker should remain conservative because the Rust code uses local helper
abstractions instead of a one-to-one inline source shape:

1. `vdw_radius(...)`
   - Source role: `PeriodicTable::getRvdw(...)` dependency
   - Verdict: keep `RDKit❗✔️`.
   - Reason: DG bounds depends on fixed local radii mapping rather than the
     original RDKit periodic-table object.

2. `ideal_bond_angle(...)` / `set_ring_angle(...)`
   - Source role: `_setRingAngle(...)` dependency
   - Verdict: keep `RDKit❗✔️`.
   - Reason: behavior is aligned for the modeled hybridization/ring-size space,
     but expressed through shared local geometry helpers.

3. `set_13_bounds_helper(...)`
   - Source role: `_set13BoundsHelper(...)`
   - Verdict: keep `RDKit❗✔️`.
   - Reason: the law-of-cosines behavior is source-backed, but the helper still
     sits behind local abstraction rather than a line-for-line low-level port.

4. `Path14Configuration`, `Path14Type`, `ComputedData`
   - Source role: local data structures mirroring RDKit state
   - Verdict: keep `RDKit❗✔️`.
   - Reason: these are representation-level mirrors, not direct RDKit function
     bodies; their behavior is only meaningful through the tested call sites.

5. `get_atom_stereo(...)`
   - Source role: `DGeomHelpers::_getAtomStereo(...)`
   - Verdict: keep `RDKit✔️✔️`.
   - Reason: the helper now carries the direct copied RDKit source body, the
     XOR-based flip rule is implemented exactly for the modeled bond stereo
     state space, and the local implementation shape is complexity-equivalent.

6. `set_lower_bound_vdw(...)`
   - Source role: `setLowerBoundVDW(...)`
   - Verdict: keep `RDKit❗✔️`.
   - Reason: observed DG bounds behavior is aligned for the audited paths, but
     the helper still relies on local donor/acceptor and radii abstractions.

7. `set_14_bounds(...)`
   - Source role: `set14Bounds(...)`
   - Verdict: keep `RDKit❗✔️`.
   - Reason: the dispatch graph is source-backed and tested, but the Rust body
     is still summarized at the top-level frame and delegates through local
     containers/helpers instead of a tighter inlined source shape.

8. `set_15_bounds(...)`
   - Source role: `set15Bounds(...)`
   - Verdict: keep `RDKit❗✔️`.
   - Reason: same as `set_14_bounds(...)`: the functional path is covered, but
     the top-level translation is still expressed through local helper/storage
     abstractions.

## Audit Outcome

For the selected DG bounds baseline, the exact closure state is:

1. no remaining known first-axis `❌` gaps in the active call chain,
2. remaining `✔️❌` markers are deliberate performance-cost disclosures,
3. remaining `✔️❗` markers are deliberate unresolved-performance disclosures,
4. remaining `❗✔️` markers are conservative abstraction/source-framing
   disclosures,
5. therefore this scope should be documented as behaviorally source-backed with
   residual marker caveats, not as blanket zero-gap parity.

## Inventory And Validation Consistency Check

The final closure statement above is allowed only because the post-Step-159
inventory and validation evidence are consistent:

1. `dev/gap_reports/dg_bounds_remaining_source_scan.md` now lists
   `DGeomHelpers::_getAtomStereo(...)` as a direct baseline function.
2. `crates/cosmolkit-core/src/chemistry/distgeom.rs` now carries copied RDKit
   source lines for `_getAtomStereo(...)`.
3. Step 157's literal checklist command was rerun and is still invalid in
   `cargo test`:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict distgeom get_atom_stereo -- --nocapture
```

Observed result:

- `error: unexpected argument 'get_atom_stereo' found`

4. The current targeted strict validation for that baseline addition passed via
   the closest executable strict command:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict get_atom_stereo -- --nocapture
```

Observed result:

- `3 passed; 0 failed`
- `1003 filtered out`
