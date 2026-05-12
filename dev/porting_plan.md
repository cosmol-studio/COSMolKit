# COSMolKit Source Porting Plan

## Navigation

| Document | Purpose |
|----------|---------|
| [porting_inventory.md](porting_inventory.md) | Per-feature code-first status (implemented / substantial / partial / deferred) |
| [work.md](work.md) | Active phase — current task, design detail, and ordered backlog |

## Overview

This plan tracks the work required for the redesigned Rust core
(`crates/cosmolkit-core/`) to reach the features that `README.md` already marks
as complete. It covers only checked README items; unchecked roadmap items are
out of scope.

Current planning mode is **code-first feature closure**. The current source has
advanced ahead of some older `dev/` status notes, so agents must treat Rust code
as the primary evidence and use Markdown only as a ledger. Strict RDKit parity
hardening and Python-facing work are deferred for now.

The plan is organized into numbered phases. Each phase has a concrete
deliverable and is executed sequentially. The current phase is tracked in
[work.md](work.md); when a phase completes, the next one from the backlog is
promoted automatically.

## Agent Entry Contract

This file is the task entry point for agents continuing the source-porting
work. `dev/README.md` is the developer operating manual; it defines rules and
commands, but the active task is defined here and in `work.md`.

Current objective:

```text
Make the redesigned Rust core honestly satisfy every feature that README.md
currently claims is already available, excluding Python bindings and strict
parity hardening for this pass.
```

First task for every new continuation:

1. Check README claims against current Rust source behavior, support specs,
   operation bodies, and then `porting_inventory.md`.
2. Fix any over-complete or stale status claim before implementing a later
   README feature. If the mismatch is caused by code pretending support exists,
   fix the code path or make it fail explicitly; Markdown edits alone are not a
   valid completion.
3. Keep `work.md` as the immediate execution queue and `porting_inventory.md`
   as the status ledger.

Do not mark a README claim complete merely because a skeleton, old-core
implementation, partial subset, or experimental branch exists. Unsupported or
partial behavior must remain visible in the inventory and must fail explicitly
in code.

Phase 0 is a code gate, not a documentation cleanup phase. An agent must not
advance past Phase 0 after only editing Markdown. It may update Markdown to
record an audit result, but Phase 0 remains active until the audited source
paths either implement the claim or reject unsupported behavior through the
operation/support contracts.

## Scope Rules

1. **Code-first feature closure.** Close functional gaps in the redesigned Rust
   core first. Do not start a strict RDKit parity campaign, add new RDKit
   comparison fixtures, or treat bit-identical parity as a blocker in this pass.
   Existing tests may still be run for regression and contract validation.

2. **Rust-only.** Python bindings, Python tests, Python stubs, packaging, and
   Python documentation are explicitly out of scope. When a checked README item
   mixes Rust and Python behavior, execute only the Rust portion now.

3. **Use pinned upstream source when filling chemistry behavior.** Keep existing
   source-line markers aligned with
   [`source_reproduction_protocol.md`](source_reproduction_protocol.md), but
   defer strict parity hardening unless the missing feature cannot be completed
   safely without it.

4. **Add typed state when needed.** If the redesigned core lacks a field needed
   by the source behavior, add it when it fits project policy: value-style API,
   no public mutable storage, operation-controlled mutation.

5. **No old-core dependency.** Do not depend on `crates/cosmolkit-core-old/`,
   re-export from it, or silently copy old-core behavior.

6. **Reread normative documents** before each implementation step. At minimum:
   - `dev/README.md`
   - `dev/policy_invariants.md`
   - `dev/topology_operations.md` for molecule operation edits
   - `dev/Macro-ControlledStateMigrationDesign.md` for operation registry edits
   - `dev/testing_contract.md` for validation policy

7. **Topology changes go through registered operations** and `OpParts`. Free
   functions may compute assignments but must not mutate `Molecule` directly.

8. **Unsupported branches must fail explicitly** with a structured error. They
   must not silently no-op or guess chemistry.

9. **Final validation** for any core or operation edit:
   ```bash
   cargo fmt --all
   cargo check -p cosmolkit-core --features op-contracts-strict
   cargo test -p cosmolkit-core --features op-contracts-strict
   ```

## Phase Map

| # | Area | Key Deliverables | Status |
|---|------|------------------|--------|
| 0 | Code Baseline Sync | Audit README claims against current Rust source; update ledger from code evidence, not old Markdown | Active |
| 1 | SMILES Writer Closure | Remove remaining unsupported writer branches needed for README-level options | Partial / substantial |
| 2 | Morgan Fingerprint Closure | Finish README-level options and document non-parity hash semantics | Partial / substantial |
| 3 | SMILES Parser Closure | Close user-visible parser gaps and sanitize integration gaps | Partial / substantial |
| 4 | MOL/SDF Writer Closure | Finish README-level writer branches and unsupported query/stereo cases | Partial / substantial |
| 5 | MOL/SDF Reader Closure | Finish V3000 and robust multi-record reader gaps | Partial |
| 6 | DG Bounds Closure | Stabilize existing bounds matrix implementation and public support status | Partial |
| 7 | 2D Coordinates Closure | Stabilize existing coordinate generation and unsupported branches | Partial / substantial |
| 8 | Drawing Closure | Stabilize existing SVG/PNG renderer, font path, and image export behavior | Partial / substantial |
| 9 | Chemistry Core Edge Closure | Hydrogens, valence/radicals, rings, aromaticity, sanitize edge branches | Partial |
| 10 | Stereochemistry Closure | CIP/geometric perception and E/Z behavior needed by current APIs | Partial |
| 11 | Batch Operations Closure | n_jobs, scheduling, transforms, masks/reports, SDF/image export | Partial |

## README Claim Closure Order

Use this order after Phase 0 has synchronized the ledger with code-level
evidence. Phase 0 is not clean merely because plan files were edited.

1. Synchronize `porting_inventory.md` and `work.md` with current code evidence.
2. Close remaining SMILES writer unsupported branches needed by README-level
   options.
3. Close Morgan fingerprint option gaps, while deferring strict hash parity.
4. Stabilize 2D coordinates, drawing/SVG/PNG, and image export.
5. Stabilize DG bounds and support status.
6. Close MOL/SDF reader and writer remaining functional branches.
7. Close hydrogen, sanitize, valence, ring, aromaticity, and stereochemistry
   functional gaps.
8. Close batch-native Rust workflows: scheduling, transforms, errors,
   masks/reports, SDF/image export.

The detailed execution queue lives in [`work.md`](work.md).

## Agent Workflow

When entering this project for porting work:

1. Read this file to understand the scope and rules.
2. Open [work.md](work.md) for the current task and design detail.
3. After completing the current phase, follow the Completion Protocol in
   `work.md` to auto-advance to the next phase.
