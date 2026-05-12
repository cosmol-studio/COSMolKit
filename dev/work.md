# Current Work

This is the immediate execution queue for the source-porting plan in
[`porting_plan.md`](porting_plan.md). The goal is to make the redesigned Rust
core honestly satisfy every feature that `README.md` currently claims is
already available.

`dev/README.md` is the developer operating manual. This file and
`porting_plan.md` are the agent task entry points.

## Current State

Use current Rust code as the primary evidence. Older Markdown status notes may
lag behind the implementation and must not be treated as authoritative.

Current code-first estimate:

- Redesigned Rust core completion against checked README claims, excluding
  Python and strict parity hardening: about 65-75%.
- Working planning number: about 70% complete, about 30% remaining.
- The remaining work is mostly functional closure, unsupported-branch removal,
  API/status consistency, and edge-case hardening. Strict RDKit parity and
  Python-facing work are deferred for now.

Recent audit results:

- `MoleculeReadParts` raw molecule access was blocked by moving it into
  `read_parts.rs` with a private field and no `molecule()` escape hatch.
- Source files with RDKit/Gemmi two-axis markers reference
  `dev/source_reproduction_protocol.md`.
- `cargo check -p cosmolkit-core --features op-contracts-strict` and
  `cargo test -p cosmolkit-core --features op-contracts-strict` passed before
  the latest documentation-only plan edits.
- Current source has moved ahead of older `porting_inventory.md` entries for
  DG bounds, drawing/SVG/PNG, Morgan, SMILES writer, SDF, and batch.

## Task 1 — Phase 0 Code Baseline Sync

**Status:** active; cannot be completed by Markdown-only edits

Phase 0 is a code gate. Its purpose is to synchronize the plan with current
Rust implementation, then immediately continue with functional code work. A new
agent must not complete Phase 0 by only editing Markdown.

Required checks:

1. For each checked README feature, inspect the current Rust entry point,
   feature support spec, operation body where applicable, unsupported branches,
   and tests.
2. Update `porting_inventory.md` only after source inspection.
3. If a feature spec or public API falsely implies support, fix the Rust code or
   make unsupported behavior fail explicitly.
4. Keep Python out of scope for this pass.
5. Do not start strict RDKit parity work in this pass. Record parity gaps as
   deferred unless they block functional correctness.

Useful audit commands:

```bash
rg -n "unsupported_stage|UnsupportedFeature|todo!|unimplemented!" \
  crates/cosmolkit-core/src -g '*.rs'

rg -n "FeatureStatus::Supported|FeatureStatus::Experimental|pub const .*_FEATURE" \
  crates/cosmolkit-core/src/support.rs
```

Done criteria:

- `porting_plan.md`, `porting_inventory.md`, and this file agree with current
  Rust code.
- Every checked README Rust/core feature has a current status and concrete
  functional gaps.
- No later task is blocked on stale Markdown saying a code path is unsupported
  when implementation exists.
- No feature is treated as complete while public code still has README-relevant
  unsupported branches.
- If Rust code changed, run:

  ```bash
  cargo fmt --all
  cargo check -p cosmolkit-core --features op-contracts-strict
  cargo test -p cosmolkit-core --features op-contracts-strict
  ```

## README Claim Backlog

The following tasks close checked README claims for the redesigned Rust core.
They intentionally exclude Python binding work and strict RDKit parity
hardening for now.

### 2. SMILES Writer Functional Closure

Current source evidence:

- `smiles_write.rs` implements plain/canonical/noncanonical output, CXSMILES,
  fragment APIs, random SMILES, atom/bond formatting, aromatic output, stereo
  and directional bond branches.
- Several `unsupported_stage(...)` branches remain.

Remaining work:

- Audit every remaining `unsupported_stage` and classify it as README-relevant
  or deferred.
- Implement README-relevant atom and bond writer branches.
- Keep full-molecule kekulization and topology changes routed through
  registered operations.

### 3. Morgan Fingerprint Functional Closure

Current source evidence:

- `fingerprint.rs` implements bit vectors, Tanimoto, Morgan output,
  count-simulation, chirality, feature invariants, custom atom/bond invariants,
  and additional output structures.
- Strict bit-identical RDKit hash parity is not the current target.

Remaining work:

- Audit option coverage against README claims.
- Ensure unsupported or intentionally deferred options fail explicitly.
- Document non-parity hash semantics in support docs if they remain by design.

### 4. 2D Coordinates And Drawing Closure

Current source evidence:

- `coordinates.rs` and operation code provide 2D coordinate generation through
  operation-controlled mutation.
- `draw.rs` implements SVG/PNG rendering, embedded Noto Sans font loading,
  atom labels, bond geometry, radical dots, scale calculation, and PNG
  rasterization through Rust libraries.
- `batch.rs` can call drawing for image export.

Remaining work:

- Audit coordinate unsupported branches and error behavior.
- Add functional checks for representative drawing outputs and PNG byte
  generation where missing.
- Ensure support specs and inventory reflect implemented drawing status.
- Defer visual parity hardening.

### 5. DG Bounds Closure

Current source evidence:

- `distgeom.rs` contains `dg_bounds_matrix` with bond/angle/torsion/VDW bounds,
  smoothing, and local tests.
- Older inventory entries that call DG unsupported are stale.

Remaining work:

- Audit public support status and error behavior.
- Add or adjust functional tests for common molecule classes if needed.
- Record strict RDKit distance-matrix parity as deferred.

### 6. MOL/SDF Reader/Writer Closure

Current source evidence:

- `io/sdf.rs` contains SDF reader, dataset indexing, MolBlock parsing, and
  V2000/V3000 writer logic with many query/stereo/property branches.
- Some unsupported branches remain for complex query predicates and reader
  edge cases.

Remaining work:

- Close V3000 reader and multi-record robustness gaps needed by README.
- Audit remaining unsupported query/stereo/SGroup branches.
- Preserve typed query, stereo, SGroup, and property state; do not flatten into
  string-only behavior.

### 7. Chemistry Core Edge Closure

Current source evidence:

- `hydrogens.rs`, `sanitize.rs`, `valence.rs`, `rings.rs`, `aromaticity.rs`,
  `kekulize.rs`, `stereo.rs`, and `ops.rs` have substantial operation-routed
  implementations and tests.

Remaining work:

- Close README-visible hydrogen, sanitize, valence, ring, aromaticity, and
  stereo gaps.
- Keep topology mutations behind registered operations and `OpParts`.
- Treat CIP/geometric stereo perception as functional work only where the
  current public API needs it; full parity hardening is deferred.

### 8. Batch Workflow Closure

Current source evidence:

- `batch.rs` supports ordered construction from SMILES, valid masks, error
  records, filtering, registered-operation transforms, SMILES export, and image
  export.
- Full `n_jobs` scheduling, richer reports, selection API breadth, and parallel
  SDF/image export still need closure.

Remaining work:

- Implement or explicitly defer `n_jobs` and scheduler behavior.
- Complete masks/reports/selection APIs needed by README-level Rust behavior.
- Complete SDF/image export behavior on top of existing writer/drawing code.

## Completion Protocol

When a task is completed:

1. Update `porting_inventory.md` from source evidence.
2. Update this file so the next task is unambiguous.
3. Run strict validation if Rust code changed:

   ```bash
   cargo fmt --all
   cargo check -p cosmolkit-core --features op-contracts-strict
   cargo test -p cosmolkit-core --features op-contracts-strict
   ```

4. Do not mark a checked README Rust/core claim complete until its public code
   path works or unsupported behavior is explicit and documented as deferred.
