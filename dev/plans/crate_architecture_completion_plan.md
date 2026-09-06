# Crate Architecture Completion Plan

## Execution Contract

- This plan is for sequential continuous execution.
- "Sequential continuous execution" means execute one step at a time in order and continue to the next unchecked step until all steps are completed, blocked, or the user interrupts.
- It does not mean executing steps in unordered batches or postponing validation for a batch of changes.
- Execute unchecked steps in order.
- Continue executing the plan until all steps are completed, blocked, or the user interrupts.
- Do not stop after every step unless the plan explicitly says to stop.
- Mark each completed step by changing only its `[ ]` to `[x]`.
- Never execute unchecked steps out of order.
- Never summarize, skip, or reinterpret later unchecked steps.
- Never treat a required reading step as “already read”.
- Do not assume the agent is diligent.
- Do not assume the model context is long enough.
- Do not rely on memory from previous turns when a required reading step is present.
- Every real task step must be immediately preceded by reading `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
- The reading step must explicitly reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
- `Implement`, `Port`, `Modify`, `Update`, and `Fix` steps must produce a concrete artifact.
- `Audit` steps must produce a written gap report and must not replace implementation steps.
- If a step adds or updates tests, the next real task after the required reading step must run the most specific relevant test command for those tests.
- Do not defer tests added for one behavior to a final whole-plan validation step.
- Final whole-plan validation is still required when the plan changes code, but it does not replace immediate targeted validation after test-writing steps.
- If the plan violates this contract, regenerate the plan before doing any work.
- Copying C++ comments, adding a dispatch stub, or adding placeholder branches is not a completed `Port` step.
- Do not use “smallest subpart”, skeleton code, dispatch-only code, placeholder code, TODO-only code, or partial porting as completion.
- If a step is too large to complete fully, regenerate the plan by splitting that step into smaller full-port steps.
- Each split step must still represent a complete source-backed behavior, not a placeholder.
- A `Port` step is complete only when the selected RDKit behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit every remaining `cosmolkit-core` dependency on `Molecule`, `MoleculeBuilder`, `OpParts`, the operation registry, contract matrices, and runtime caches and write `dev/gap_reports/crate_architecture_completion_inventory.md` with owner, callers, feature gate, source-port status, and migration order for each item.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Update `dev/crate_architecture.md` with the accepted production transaction shape proven by `dev/experiments/crate_split`, including direct `OpParts<'_, Access>` and `MultiOutputOpParts<'_, Access>` operation bodies and direct multi-output block tuples.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Audit `cosmolkit-types` and `cosmolkit-model` against all existing production value-level invariants and add a gap table for unresolved ID, adjacency, topology, coordinate, and property alignment checks to `dev/gap_reports/crate_architecture_completion_inventory.md`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Implement the first complete missing model-level invariant family identified in the inventory without introducing a runtime dependency.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Add focused model-crate tests for the invariant family implemented in Step 9.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Run the focused `cosmolkit-model` test command for the tests added in Step 11.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Complete the remaining model-level invariant families from the inventory one complete family at a time and record the completed coverage in the inventory.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Add or update focused `cosmolkit-model` tests for every invariant family completed in Step 15.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Run the focused `cosmolkit-model` test commands for the tests added in Step 17.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Audit the valence, sanitization, and hydrogen implementation boundary and write the complete migration scope, source dependencies, and required detached-block interfaces to `dev/gap_reports/crate_architecture_completion_inventory.md`.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Implement the complete read-only valence assignment slice behind an explicit detached-block API in `cosmolkit-core`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Add source-backed parity tests for the detached valence assignment slice implemented in Step 23.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Run the focused strict core parity tests added in Step 25.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Audit the runtime ownership migration and write the exact module move order, dependency cut points, and generated-code constraints to `dev/gap_reports/crate_architecture_completion_inventory.md`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Audit the complete hydrogen planner and runtime-apply boundary and write its detached inputs, result records, and remaining runtime dependencies to `dev/gap_reports/crate_architecture_completion_inventory.md`.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Implement the complete detached AddHs planner using explicit model blocks and a value valence assignment.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Add source-backed parity tests for the detached AddHs planner.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run the focused strict AddHs parity tests added in Step 35.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Audit the complete RemoveHs planner into detached planning and runtime application boundaries and write the split to `dev/gap_reports/crate_architecture_completion_inventory.md`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Implement the complete baseline RemoveHs candidate-selection slice over explicit detached topology and properties.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add source-backed parity tests for the baseline detached RemoveHs candidate-selection slice.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run the focused strict baseline RemoveHs parity tests added in Step 43.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Implement the complete detached RemoveHs stereo, isotope, SGroup, and sanitize planning slices over explicit model values.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Add source-backed parity tests for the detached RemoveHs planning slices added in Step 47.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run the strict RemoveHs planning parity tests added in Step 49.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Implement the complete detached property-cache sanitization slice while leaving cache installation and operation lifecycle in the runtime.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Add source-backed parity tests for the detached property-cache sanitization slice.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run the focused strict sanitization parity tests added in Step 55.
Step 58 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [ ]: Move search, SMARTS, CX, notation, and IO algorithms behind complete explicit-model APIs with no live runtime imports in their implementation paths.
Step 60 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [ ]: Add source-backed parity tests for the search, notation, and IO APIs migrated in Step 59.
Step 62 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [ ]: Run the strict search, notation, and IO parity tests added in Step 61.
Step 64 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [ ]: Move descriptors and fingerprints behind complete explicit-model APIs with no live runtime imports in their implementation paths.
Step 66 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [ ]: Add source-backed parity tests for descriptors and fingerprints migrated in Step 65.
Step 68 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [ ]: Run the strict descriptor and fingerprint parity tests added in Step 67.
Step 70 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [ ]: Move coupled 3D, force-field, alignment, stereoisomer, and tautomer algorithms behind complete explicit-model APIs with no live runtime imports in their implementation paths.
Step 72 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [ ]: Add source-backed parity tests for the coupled algorithm families migrated in Step 71.
Step 74 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [ ]: Run the strict parity tests added in Step 73.
Step 76 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [ ]: Move the complete `Molecule` and `MoleculeBuilder` implementation into `cosmolkit` while preserving one authoritative runtime and compiling all current operation callers.
Step 78 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [ ]: Move the complete `OpParts`, operation registry, generated matrices, runtime cache authority, and contract lifecycle into `cosmolkit`.
Step 80 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [ ]: Update every generated operation wrapper and body to use the `cosmolkit` runtime’s access-marker parameterized `OpParts` or `MultiOutputOpParts` directly.
Step 82 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [ ]: Add strict runtime tests covering detached single-output rollback, in-place rollback, capability-restricted block visibility, multi-output candidate validation, and final commit validation.
Step 84 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [ ]: Run the focused strict runtime tests added in Step 83.
Step 86 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [ ]: Rebuild Cargo fine-grained feature gates and user-facing bundles so every public capability selects only its own domain dependency and no implementation edge exposes an unrelated `Molecule` method.
Step 88 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [ ]: Add feature-matrix compile tests proving individual features and bundles expose exactly their declared public APIs.
Step 90 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [ ]: Run the feature-matrix compile tests added in Step 89.
Step 92 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [ ]: Update Rust examples, Python bindings, documentation, integration tests, and package metadata so supported users import `Molecule` only from `cosmolkit` and no public path relies on `cosmolkit_core::Molecule`.
Step 94 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [ ]: Add cross-language compatibility tests for the public paths updated in Step 93.
Step 96 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [ ]: Run the most specific Rust and Python compatibility tests added in Step 95.
Step 98 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [ ]: Remove transitional core runtime exports, duplicate adapters, double-underscore migration feature aliases made obsolete by split crates, and every remaining public `cosmolkit_core::Molecule` path.
Step 100 [ ]: Write `dev/gap_reports/crate_architecture_completion_report.md` proving the final dependency graph is acyclic below `cosmolkit`, `Molecule` and the sole runtime are owned by `cosmolkit`, and no child crate has runtime or derived-cache authority.
Step 101 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 102 [ ]: Run `cargo fmt --all`, `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`, the feature-matrix checks, and the Python binding validation required by `AGENTS.md`.

## Current execution checkpoint

Step 59 is in progress. The search predicate/context slice now uses the
detached `SearchTargetAccess` boundary and has a block-only context constructor
with a strict focused test. Live `Molecule` matcher entrypoints remain
transitional adapters; Step 59 is not complete until matcher, SMARTS writer,
CX/notation, and IO implementation paths have the same boundary and those
adapters are owned by `cosmolkit`.
