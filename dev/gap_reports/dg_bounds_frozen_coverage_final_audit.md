# DG Bounds Frozen Coverage Final Audit

Date: 2026-05-18

## Final Conclusion

For the DG bounds baseline defined by
`dev/dg_bounds_rdkit_full_port_checklist.md`, coverage is exactly `100%`.

That `100%` statement is scoped and means:

1. every baseline direct RDKit DG bounds helper/function has a source-backed
   Rust implementation,
2. the direct baseline inventory now explicitly includes
   `DGeomHelpers::_getAtomStereo(...)` and its copied source anchor,
3. the current targeted strict validation for that added direct helper passes,
4. the previously blocking matrix-layout, triangle-smoothing, and wrapper
   default gaps are closed,
5. the local strict DG bounds validation slice passes.

It does **not** mean the entire `cosmolkit-core` crate is globally free of
non-DG failures.

## Audited Source Baseline

- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`
- `third_party/rdkit/Code/DistGeom/BoundsMatrix.h`
- `third_party/rdkit/Code/DistGeom/TriangleSmooth.cpp`
- `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp`

## Rust Artifacts Audited

- `crates/cosmolkit-core/src/chemistry/distgeom.rs`
- `dev/gap_reports/dg_bounds_remaining_source_scan.md`
- `dev/porting_inventory.md`
- `crates/cosmolkit-core/src/support.rs`

## Frozen Coverage Verdict

The frozen DG bounds verdict is:

- baseline direct helper coverage: complete
- source reproduction closure: complete within the selected DG bounds scope
- local DG strict validation: passing
- exact DG bounds coverage claim: `100%`

## Evidence

### Source rescan

`dev/gap_reports/dg_bounds_remaining_source_scan.md` now records zero remaining
direct first-axis gaps for the selected DG bounds baseline, including the
explicit `_getAtomStereo(...)` entry added in the final baseline refresh.

### Local DG validation

Executed command:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict get_atom_stereo -- --nocapture
```

Result at freeze time:

- `3 passed; 0 failed`
- `1003 filtered out`

For completeness, the literal Step 157 checklist command was rerun and remains
invalid in `cargo test`:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict distgeom get_atom_stereo -- --nocapture
```

Observed result:

- `error: unexpected argument 'get_atom_stereo' found`

### Final project-level validation commands executed

Executed:

```bash
cargo fmt --all
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict
```

Observed results:

1. `cargo fmt --all`: passed
2. `cargo check -p cosmolkit-core --features op-contracts-strict`: passed
3. `cargo test -p cosmolkit-core --features op-contracts-strict`: failed, but
   the remaining failures are outside the DG bounds baseline

## Non-DG Global Test Failures At Freeze Time

The full strict crate test run still reports unrelated failures in:

1. `chemistry::stereo_enumerate::tests::test_find_tetrahedral_centers`
2. `chemistry::stereo_enumerate::tests::test_enum_stereoisomers_basic`
3. `chemistry::stereo_enumerate::tests::test_count_functions`
4. `notation::smiles::tests::close_mol_rings_uses_second_partial_bond_cx_smiles_bond_idx_like_rdkit`
5. `notation::smiles::tests::from_smiles_assigns_ring_closure_double_bond_stereo_atoms_like_rdkit_row_86`
6. `notation::smiles::tests::from_smiles_with_sanitize_false_parses_percent_ring_numbers`
7. `notation::smiles::tests::from_smiles_with_sanitize_false_uses_closing_explicit_ring_bond_when_opening_unspecified`
8. `notation::smiles_write::tests::all_primary_smiles_writer_modes_fail_closed_until_ported`
9. `notation::smiles_write::tests::cx_no_cx_data_returns_plain_smiles`
10. `notation::smiles_write::tests::plain_smiles_writer_strips_dative_bonds_on_working_copy_when_requested`
11. `notation::smiles_write::tests::writer_canonicalizes_fused_ring_double_bond_directions_like_rdkit`

These failures block a whole-crate green build, but they do not reopen a DG
bounds source gap according to the selected checklist baseline.

## Scope Boundary Notes

This frozen `100%` claim intentionally does not expand beyond the audited DG
bounds baseline. In particular, it does not claim closure for:

1. unrelated stereo-enumeration functionality
2. unrelated SMILES parser behavior
3. unrelated SMILES writer behavior
4. broader crate-wide RDKit parity outside the DG bounds scope

If the DG bounds baseline is later expanded to include additional convenience
APIs or new RDKit surfaces, this frozen audit must be regenerated.
