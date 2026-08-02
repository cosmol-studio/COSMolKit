# Source-Port Duplicate Consolidation Plan

## Execution Contract

- This plan is for sequential continuous execution.
- Execute unchecked steps in order and continue until all steps are completed, blocked, or the user interrupts.
- Mark a completed step by changing only its `[ ]` to `[x]`.
- Every real task is immediately preceded by a fresh reading of `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
- Consolidate only behavior proven to reproduce the same upstream source function and active source branch.
- Preserve separate adapters, errors, ownership boundaries, and state mutation when those correspond to different upstream behavior.
- Do not consolidate merely because two Rust bodies look alike.
- Do not execute Git operations.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Update `dev/gap_reports/source_port_duplicate_implementation_audit.md` with an explicit same-source consolidation matrix and frozen different-source candidates.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Consolidate the four `countSwapsToInterconvert` ports into one source-backed generic kernel while preserving caller-specific failure behavior.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add focused tests for the shared `countSwapsToInterconvert` kernel and all caller failure mappings.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run the most specific `countSwapsToInterconvert` focused tests with strict operation contracts.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Consolidate the two `DGeomHelpers::_setRingAngle` ports into one source-backed pure kernel while retaining caller-specific molecule access.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Add focused branch and boundary tests for the shared `_setRingAngle` kernel.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run the most specific `_setRingAngle` focused tests with strict operation contracts.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Replace format-local copies of the same RDKit periodic-table symbol and atomic-number functions with the canonical source-backed conversion functions while preserving format errors.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Add focused exhaustive periodic-table delegation tests for MOL2, PDB, and SMILES boundaries.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run the most specific periodic-table delegation tests with strict operation contracts.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Consolidate the MMFF and UFF copies of `ForceFieldsHelper::OptimizeMoleculeConfs`, its active non-threaded thread-count branch, and `OptimizeMoleculeConfsST` into one source-backed internal driver.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Add focused MMFF and UFF tests proving shared-driver result, mutation, missing-parameter, thread-count, and iteration-limit behavior.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Run the most specific MMFF and UFF shared-driver tests with strict operation contracts.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Consolidate `RingUtils::makeRingNeighborMap` and `RingUtils::pickFusedRings` into source-backed generic kernels while retaining aromaticity and kekulization selection policies.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Add focused ordering and overlap tests for both callers of the shared RingUtils kernels.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run the most specific RingUtils focused tests with strict operation contracts.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Consolidate the source-identical pure portions of `Canon::chiralAtomNeedsTagInversion`, `atomHasFourthValence`, `hasSingleHQuery`, and `Chirality::getChiralPermutation` behind explicit builder and molecule read adapters.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Add focused parser and writer tests for every branch of the consolidated Canon helpers.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run the most specific Canon-helper focused tests with strict operation contracts.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Consolidate the two `getIsoMap` ports into one source-backed traversal kernel while preserving builder and molecule index/error adapters.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add focused tests for all source branches and caller failure mappings of the shared `getIsoMap` kernel.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run the most specific `getIsoMap` focused tests with strict operation contracts.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Update the duplicate-implementation audit with completed consolidations, intentionally separate different-source paths, and remaining blockers.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run `cargo fmt --all -- --check`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.
