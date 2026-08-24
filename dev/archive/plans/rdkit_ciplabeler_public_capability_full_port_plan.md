# RDKit CIPLabeler Public Capability Full Port Plan

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
- Every real task step must be immediately preceded by reading: `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
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

## Scope And Ownership

- The source pin is RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8c`, as recorded in `testdata/reference/rdkit.json`.
- `crates/cosmolkit-core/src/chemistry/ciplabeler.rs` is the sole modern CIPLabeler implementation core to audit and complete.
- Legacy `_CIPRank` and legacy `_CIPCode` paths remain distinct and may survive only for demonstrated pinned legacy consumers.
- Modern assignment must never delegate to the legacy `Chirality.cpp` rank/code implementation.
- Shared source dependencies such as atropisomer carrier extraction must have one internal implementation owner.
- `findPotentialStereo`, full `assignStereochemistry`, `assignStereochemistryFrom3D`, and 3D stereo perception remain separate capabilities unless a source dependency is proved and recorded.
- The public support claim is limited to configurations constructed by the pinned `findConfigs()` dispatcher.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit every pinned CIPLabeler source function against `crates/cosmolkit-core/src/chemistry/ciplabeler.rs` and update `dev/gap_reports/rdkit_ciplabeler_public_capability_source_inventory.md` with a reviewed per-function status and every discovered gap.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit every modern CIPLabeler consumer and every legacy CIP rank/code entry point and record the single-owner dependency graph plus each retained or narrowed legacy public boundary in `dev/gap_reports/rdkit_ciplabeler_public_capability_source_inventory.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Audit `_CIPCode`, `_CIPNeighborOrder`, `_CIPComputed`, preexisting `_CIPRank`, bond stereo, and stereo-atom lifecycle against pinned RDKit and record the full/selected/success/failure/clone/serialization/mutation matrix in `dev/gap_reports/rdkit_ciplabeler_public_capability_source_inventory.md`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Audit pinned exception, recursion-limit, index-selection, malformed-state, and cancellation behavior and record exact Rust/Python error mappings or explicit unsupported boundaries in `dev/gap_reports/rdkit_ciplabeler_public_capability_source_inventory.md`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Define the focused oracle schema, fixture ownership, maintained corpus tiers, and complete observable-state comparison fields in `dev/gap_reports/rdkit_ciplabeler_public_capability_source_inventory.md`.

Step 12 [x]: Read `dev/README.md`.
Step 13 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 14 [x]: Reproduce the pinned `Descriptor`, `Priority`, `Sort`, and `Edge` functions in the single modern core with exact source anchors and focused unit evidence.
Step 15 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 16 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_sort_priority`.

Step 17 [x]: Read `dev/README.md`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Reproduce the pinned `Node` constructors, accessors, child creation, flags, edge access, and auxiliary-state functions in the single modern core with exact source anchors and focused unit evidence.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_node`.

Step 22 [x]: Read `dev/README.md`.
Step 23 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 24 [x]: Reproduce the pinned `Digraph` construction, lazy expansion, node lookup, Rule6 reference, and root-change functions in the single modern core with exact source anchors and focused unit evidence.
Step 25 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 26 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_digraph`.

Step 27 [x]: Read `dev/README.md`.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Reproduce the pinned `CIPMol` molecule view and all Mancude fractional-atomic-number functions in the single modern core with exact source anchors and focused unit evidence.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_cipmol`.

Step 32 [x]: Read `dev/README.md`.
Step 33 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 34 [x]: Reproduce the pinned `SequenceRule` recursion, comparison, sorting, upward-edge, and recursion-budget functions plus `Rules` aggregation in the single modern core with exact source anchors and focused unit evidence.
Step 35 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 36 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_sequence_rule`.

Step 37 [x]: Read `dev/README.md`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Reproduce pinned `Rule1a::compare`, `Rule1b::compare`, `Rule2::compare`, and `Rule3::ord/compare` in the single modern core with exact source anchors and focused unit evidence.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_rules_rule`.

Step 42 [x]: Read `dev/README.md`.
Step 43 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 44 [x]: Reproduce pinned `PairList`, `Rule4a`, `Rule4b`, and `Rule4c` functions in the single modern core with exact source anchors and focused unit evidence.
Step 45 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_pairlist`.

Step 47 [x]: Read `dev/README.md`.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Reproduce pinned `Rule5New` and `Rule6` functions in the single modern core and document old `Rule5` as unreachable from the modern assignment dispatcher with exact source anchors and focused unit evidence.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_rules_rule5new_rule6`.

Step 52 [x]: Read `dev/README.md`.
Step 53 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 54 [x]: Reproduce pinned `Configuration` base functions and `Tetrahedral` constructors, labels, parity, carrier ordering, and primary-state functions in the single modern core with exact source anchors and focused unit evidence.
Step 55 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 56 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_tetrahedral`.

Step 57 [x]: Read `dev/README.md`.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Reproduce pinned `Sp2Bond` functions and the exact reachable `Chirality::findStereoAtoms` dependency in the single modern core with exact source anchors and focused unit evidence.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_sp2bond`.

Step 62 [x]: Read `dev/README.md`.
Step 63 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 64 [x]: Consolidate pinned `Atropisomers::getAtropisomerAtomsAndBonds` carrier extraction into one exact internal owner used by both stereochemistry and modern CIP labeling.
Step 65 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 66 [x]: Add focused unit tests proving the consolidated atropisomer carrier helper reproduces endpoint and neighbor-bond ordering without a duplicate CIP-local implementation.
Step 67 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 68 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict atropisomer_atoms_and_bonds`.

Step 69 [x]: Read `dev/README.md`.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Reproduce pinned `AtropisomerBond` constructors, labels, carrier ordering, parity, and primary-state functions in the single modern core with exact source anchors and focused unit evidence.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_atropisomerbond`.

Step 74 [x]: Read `dev/README.md`.
Step 75 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 76 [x]: Reproduce pinned `findConfigs`, `labelAux`, `label`, both `assignCIPLabels` overloads, selection masks, two-pass recursion budget, and exact completion-state behavior in the single modern core with exact source anchors and focused unit evidence.
Step 77 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 78 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_assign`.

Step 79 [x]: Read `dev/README.md`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Fix all exact-width source types and sentinels in the modern core, including `Atom::NOATOM`, `_CIPNeighborOrder`, selection indexes, recursion counters, and PairList pairing storage.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Add focused boundary tests for every exact-width conversion, overflow, sentinel serialization, invalid selection index, and recursion-counter edge fixed in Step 81.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_exact_width`.

Step 86 [x]: Read `dev/README.md`.
Step 87 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 88 [x]: Implement the audited cancellation boundary through an existing real project cancellation mechanism or a structured unsupported error without claiming interrupted-call parity.
Step 89 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 90 [x]: Add focused Rust tests for the cancellation or unsupported-cancellation behavior implemented in Step 88.
Step 91 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 92 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict ciplabeler_cancellation`.

Step 93 [x]: Read `dev/README.md`.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Add the project-native typed `CipDescriptor` conversion and atom/bond descriptor query surface backed only by modern CIP property state.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Add focused public Rust tests for typed R/S/r/s, E/Z, M/P/m/p, absent descriptors, malformed stored descriptors, and the unsupported non-tetrahedral assignment boundary.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_ciplabeler_public_api`.

Step 100 [x]: Read `dev/README.md`.
Step 101 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 102 [x]: Add a `molecule_ops!` entry providing value-style and trailing-underscore in-place modern CIP assignment for full, selected-atom, selected-bond, and recursion-limit options.
Step 103 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 104 [x]: Add strict operation-contract tests for read/write capabilities, full and selected assignment, value preservation, in-place partial-failure validity, mapping, and invariant behavior.
Step 105 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 106 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict cip_label`.

Step 107 [x]: Read `dev/README.md`.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Implement the audited modern CIP state invalidation, preservation, remapping, clone, and supported serialization lifecycle across all registered topology/stereo operations.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Add a generated operation-matrix test proving the lifecycle from Step 109 for every registered topology/stereo mutation rather than a hand-selected subset.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict cip_state_lifecycle`.

Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Add the pinned-RDKit focused CIP oracle generator and register it in `tools/testdata/rdkit/generate_all.py`.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Run `uv run python tools/testdata/rdkit/generate_all.py --only ciplabeler`.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Add committed focused CIP fixtures and expected output covering upstream regressions, source branches, complete observable state, selected assignment, repeated calls, recursion limits, dummy atoms, phosphine/arsine centers, pseudoasymmetry, non-kekulizable fragments, Sp2, atropisomers, and unsupported dispatcher boundaries.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Run `uv run python tools/testdata/rdkit/generate_all.py --only ciplabeler --check`.

Step 122 [x]: Read `dev/README.md`.
Step 123 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 124 [x]: Add a public Rust integration oracle that compares every focused fixture field against the single modern CIP assignment operation.
Step 125 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 126 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_ciplabeler_parity`.

Step 127 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 128 [x]: Extend the maintained CIP oracle data generation to `smiles_small` and every active complete 5,000-row corpus profile with deterministic sharding and complete state output.
Step 129 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 130 [x]: Run the maintained CIP generator for `smiles_small` and the complete 5,000-row corpus.
Step 131 [x]: Read `dev/README.md`.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Add maintained Rust CIP parity gates for `smiles_small` and every active complete 5,000-row corpus profile.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_ciplabeler_parity`.

Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Add a reproducible modern CIP assignment phase to `dev/tools/chembl_parity/` with deterministic partitioning, full and selected branches, recursion profiles, complete state comparison, and aggregate reporting.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Run the ChEMBL 37 CIP parity phase across the complete maintained corpus with the project parallel runner.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Record the complete ChEMBL 37 CIP branch counts, comparison counts, error distribution, and mismatch result in the maintained parity report artifact.

Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Add the project-native Python modern CIP assignment and atom/bond descriptor query surface backed by the Rust operation and typed state.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Add Python tests for full assignment, atom-only selection, bond-only selection, simultaneous selection, recursion-limit errors, descriptor queries, source-state preservation, and unsupported cancellation behavior.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml`.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Run `.venv/bin/pytest python/tests -k cip`.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Regenerate `python/cosmolkit.pyi` with the project stub generator.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Run `.venv/bin/basedpyright python/tests python/examples`.

Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Add Rust and Python modern CIP assignment examples demonstrating full assignment, selected atom assignment, and typed descriptor queries without invoking legacy rank helpers.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Run all added Rust CIP examples and Python CIP examples against the built bindings.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Update current Rust API docs, Python docs, operation-system docs, support matrices, and `dev/parity_scope.md` with the exact modern CIP support boundary, lifecycle, error behavior, and measured evidence without modifying historical plans.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Run `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html`.

Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Run `cargo fmt --all -- --check`.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Run `.venv/bin/pytest`.

Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Audit all modern CIP source markers, unreproduced branches, duplicate core logic, legacy reachability, public state access, operation registry entries, fixture manifests, corpus gates, generated stubs, and documentation claims and write the final closure result into `dev/gap_reports/rdkit_ciplabeler_public_capability_source_inventory.md`.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Update this plan only by changing completed step checkboxes after every completion and move it to `dev/archive/plans/` only after no unchecked execution step remains.
