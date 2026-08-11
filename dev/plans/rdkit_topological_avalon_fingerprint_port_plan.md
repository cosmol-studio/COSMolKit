# RDKFingerprint and Avalon Source-Port Plan

## Execution Contract

- This plan is for sequential continuous execution.
- Execute unchecked steps in order and continue until all steps are completed,
  blocked, or the user interrupts.
- Mark a completed step by changing only its `[ ]` to `[x]`.
- Never summarize, skip, reorder, or reinterpret unchecked steps.
- Every real task step must be immediately preceded by reading
  `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
- Every audit must produce a concrete report, and every implementation or test
  step must produce a concrete source or test artifact.
- Do not use heuristic fingerprints, partial bitsets, similarity correlation,
  or another fingerprint family as a fallback.
- Do not upgrade support status until the complete exposed behavior is
  source-backed and exact-bit parity passes.
- During porting, do not execute any Git operation unless the user explicitly
  requests it.

## Scope and Completion Criteria

The current public methods fail closed:

- `Molecule::topological_fingerprint` / Python
  `Molecule.topological_fingerprint()`
- `Molecule::avalon_fingerprint` / Python `Molecule.avalon_fingerprint()`

The plan is complete only when both methods reproduce their pinned upstream
implementations for every exposed option, retain structured errors for any
still-unmodeled input state, and pass exact-bit focused and 5000-row parity.
Source order, integer-width behavior, random-bit generation, folding, output
size, and option interactions are behavioral requirements.

## Plan

Step 1 [ ]: Read `dev/agent_plan_standard.md`.

Step 2 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 3 [ ]: Audit the pinned RDKit RDKFingerprint wrapper, RDKitFP generator, shared fingerprint-generator helpers, path enumeration, random-bit generation, density folding, invariants, and additional-output call graph and write `dev/gap_reports/rdkit_topological_fingerprint_source_inventory.md` with exact source ranges and exposed-option mappings.

Step 4 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 5 [ ]: Implement the complete source-exact RDKit branched-path enumeration, path hashing, source random-bit generation, and density-folding behavior identified by the inventory with verbatim line-by-line two-axis source anchors.

Step 6 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 7 [ ]: Add focused exact-bit tests for RDKit path length, branched paths, bond-order selection, bits-per-hash, folding thresholds, atom inclusion, atom exclusion, disconnected fragments, rings, stereochemistry, and fixed random-bit ordering.

Step 8 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 9 [ ]: Run the focused RDKFingerprint source-parity tests with `op-contracts-strict` and require exact expected bit indexes for every row.

Step 10 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 11 [ ]: Implement the complete exposed RDKFingerprint wrapper and generator option semantics identified by the inventory without changing the public value-style API or substituting another fingerprint generator.

Step 12 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 13 [ ]: Add Rust and Python tests that cover every exposed topological-fingerprint parameter individually and in interactions, including invalid sizes and structured unsupported input states.

Step 14 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 15 [ ]: Run the focused Rust and Python topological-fingerprint API tests and require exact bit-vector and error parity.

Step 16 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 17 [ ]: Validate topological fingerprints on the maintained 5000-row corpus across the full exposed parameter matrix and record exact-bit results in `dev/gap_reports/rdkit_topological_fingerprint_full_port_validation.md`.

Step 18 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 19 [ ]: Audit the pinned Avalon/reaccs source, wrapper boundary, licensing, `bitFlags`, `isQuery`, `resetVect`, hydrogen, tautomer, atom-selection, and byte-rounded output call graph and write `dev/gap_reports/avalon_fingerprint_source_inventory.md` with exact source ranges and exposed-option mappings.

Step 20 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 21 [ ]: Implement the complete source-exact Avalon/reaccs fingerprint path identified by the inventory with verbatim line-by-line two-axis source anchors and no locally invented hash behavior.

Step 22 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 23 [ ]: Add focused exact-bit Avalon tests for every exposed option and source branch, including bit flags, query mode, reset behavior, hydrogen handling, tautomeric mode, atom selection, byte rounding, invalid sizes, and structured unsupported input states.

Step 24 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 25 [ ]: Run the focused Rust and Python Avalon API tests and require exact bit-vector and error parity.

Step 26 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 27 [ ]: Validate Avalon fingerprints on the maintained 5000-row corpus across the full exposed parameter matrix and record exact-bit results in `dev/gap_reports/avalon_fingerprint_full_port_validation.md`.

Step 28 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 29 [ ]: Update the feature registry, support matrices, README, Rust API docs, Python docs, stubs, and examples to mark only the fully validated fingerprint families as supported.

Step 30 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 31 [ ]: Run formatting, strict workspace checks, full release workspace tests, Python tests, rustdoc, Sphinx warnings-as-errors, all fingerprint examples, and the lightweight feature-combination audit as final validation.
