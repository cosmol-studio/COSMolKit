# SMILES RDKit Frozen Coverage Final Audit

Date: 2026-05-15

## Verdict

Frozen-scope coverage is **not** exactly `100%`.

The final strict validation now passes, but the frozen baseline still contains
unresolved copied-source markers in the parser, kekulize, valence, and
sanitize/orchestration paths. Under the checklist's own definition, that blocks
any `100% coverage` claim.

## Basis For Verdict

### 1. Final strict validation passed

The required validation commands for this audit segment succeeded:

- `cargo check -p cosmolkit-core --features op-contracts-strict`
- `cargo test -p cosmolkit-core --features op-contracts-strict`

The final strict test result was:

- 756 unit tests passed
- 1 doctest passed
- 0 failures

Passing strict validation is necessary for closure, but it is not sufficient on
its own because the checklist also requires a zero-marker frozen-scope rescan.

### 2. Frozen-scope marker audit is still open

The final rescan of the frozen baseline found these remaining unresolved marker
classes:

- `crates/cosmolkit-core/src/notation/smiles.rs`:
  1 `RDKit❌❌`, 2 `RDKit❗❗`, 14 `RDKit✔️❌`, 713 `RDKit❗✔️`
- `crates/cosmolkit-core/src/chemistry/kekulize.rs`:
  397 `RDKit✔️❌`
- `crates/cosmolkit-core/src/chemistry/valence.rs`:
  4 `RDKit✔️❌`
- `crates/cosmolkit-core/src/operations/ops.rs`:
  216 `RDKit✔️❌`

These frozen-scope files are now marker-closed and therefore are no longer
blocking the checklist result:

- `crates/cosmolkit-core/src/notation/smiles_write.rs`
- `crates/cosmolkit-core/src/notation/canon_rank.rs`
- `crates/cosmolkit-core/src/chemistry/rings.rs`

Because unresolved frozen-scope markers still remain elsewhere, the result
cannot be certified as `100%`.

### 3. Remaining blockers are function-level, not merely bookkeeping

The current Step 219 source scan already records the exact remaining blocking
functions. The principal open groups are:

- parser/helper functions in `notation/smiles.rs`, including `next_token`,
  `parse_mol`, `mol_from_smiles`, `cleanup_after_parsing`, ring-closure
  helpers, and several stereo/query helpers
- all frozen kekulize helpers in `chemistry/kekulize.rs`
- `ValenceContext::new` in `chemistry/valence.rs`
- sanitize/property/cleanup orchestration helpers in `operations/ops.rs`,
  including `run_sanitize_pipeline` and the cleanup assignment helpers

Until those markers are closed, final audit language must remain conservative.

## Files Updated In This Final Audit Segment

- `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`
- `dev/gap_reports/smiles_rdkit_frozen_coverage_final_audit.md`
- `dev/porting_inventory.md`
- `crates/cosmolkit-core/src/support.rs`

## Final Statement

The final strict validation pass succeeded, and several previously open
frozen-scope files are now marker-closed. Even so, frozen-scope coverage is
still **below** `100%` because unresolved parser, kekulize, valence, and
sanitize/orchestration markers remain in the frozen baseline.
