# RDKFingerprint And Avalon Full Source-Port Plan

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
- A `Port` step is complete only when the selected upstream behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Source Pins

- Chemistry oracle and RDKit source: RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- RDKFingerprint source boundary: `RDKFingerprintMol()` and its active `RDKitFPGenerator`, `FingerprintGenerator`, subgraph enumeration, hashing, random-bit expansion, folding, and `AdditionalOutput` call graph.
- Avalon adapter source boundary: the pinned RDKit `External/AvalonTools` wrapper from the same RDKit revision.
- Avalon engine source: `rdkit/ava-formake` tag `AvalonToolkit_2.0.5-pre.3`, archive MD5 `7a20c25a7e79f3344e0f9f49afa03351`.
- The Avalon archive is materialized only under ignored `tmp/parity-audit/sources/`; source provenance and licensing are recorded in a committed gap report.

## Scope

The two currently fail-closed public families are:

- `Molecule::topological_fingerprint()` and Python `Molecule.topological_fingerprint()` for RDKit `RDKFingerprintMol` behavior;
- `Molecule::avalon_fingerprint()` and Python `Molecule.avalon_fingerprint()` for RDKit `pyAvalonTools.GetAvalonFP(ROMol, ...)` behavior.

The current placeholder parameter types are not accepted as the final API contract. The RDKFingerprint placeholder omits source options and contains `ignore_atoms`, which is not an `RDKFingerprintMol` parameter. The Avalon placeholder incorrectly exposes path-fingerprint options instead of the source `nBits`, `isQuery`, and `bitFlags` surface. These shapes must be corrected while the methods still fail closed, before implementation is presented as supported.

The public API remains project-native and value-style:

- RDKFingerprint returns a fresh explicit bit vector, with a separate typed output when source `atomBits` or `bitInfo` provenance is requested;
- Avalon returns a fresh explicit bit vector, so source `resetVect` is tested inside the adapter but is not exposed as a meaningless user option for a newly allocated result;
- Avalon Python defaults follow the pinned RDKit Python wrapper, including its default bit-flag profile, rather than silently inheriting the C++ overload default;
- unsupported query or molecule states return structured errors and never approximate bits.

This plan does not port RDKit LayeredFingerprint, PatternFingerprint, unfolded/count RDK fingerprints, the public fingerprint-generator class hierarchy, Avalon count fingerprints, Avalon string-input overloads, CheckMol, canonical SMILES, or Avalon coordinate generation. Helpers reached by the two selected bit-vector call graphs remain in scope even when they live in shared upstream files.

## Completion Criteria

The plan is complete only when:

- every copied upstream line in the active call graphs has an individual two-axis source marker;
- integer widths, signedness, overflow, source iteration order, Boost hashing, Boost random-number generation, folding, byte packing, and option interactions reproduce the pinned sources;
- RDKFingerprint exact bits and requested provenance match on focused branch fixtures and every active 5,000-row profile;
- Avalon exact bits match for every public bit-flag profile, query mode, output-size boundary, aromaticity pass, and every active 5,000-row profile;
- no aggregate similarity, correlation, approximate vector, fallback fingerprint family, filtered mismatch, or one-row smoke test is accepted as parity;
- each family has a zero-unexplained-mismatch validation report with reproducible commands, timings, source versions, corpus checksums, and machine-readable artifacts under ignored `tmp/parity-audit/`;
- support metadata is family-specific and upgraded only after that family's report proves completion;
- Rust, Python, documentation, examples, stubs, feature combinations, strict tests, and release tests agree on the final surface;
- deterministic benchmarks show no unexplained algorithmic-complexity regression against the pinned source implementation.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.

Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 3 [x]: Audit the pinned `RDKFingerprintMol` active call graph, current COSMolKit placeholder surface, reusable Morgan generator machinery, reference tests, option defaults, additional outputs, and error branches and write `dev/gap_reports/rdkit_topological_fingerprint_source_inventory.md` with exact source ranges, a line-coverage ledger, and a public-option mapping.

Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 5 [x]: Modify the still-unsupported Rust and Python RDKFingerprint parameter and result types to match the selected source boundary, remove invented legacy-wrapper options, add missing source options, and retain structured fail-closed behavior.

Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 7 [x]: Add focused Rust tests for the corrected RDKFingerprint parameter defaults, invalid ranges, typed provenance allocation, and structured unsupported result.

Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 9 [x]: Run the exact focused Rust test target added in Step 7 with `op-contracts-strict`.

Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 11 [x]: Port the complete active RDKitFP atom-invariant and bond-hash input preparation source unit identified in Step 3 with verbatim line-by-line two-axis source anchors.

Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 13 [x]: Add focused exact-value tests for RDKitFP atom invariants, aromatic flags, query-bond classification, bond-order modes, custom invariants, and C++ integer-width boundaries.

Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 15 [x]: Run the exact focused Rust test target added in Step 13 with `op-contracts-strict`.

Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 17 [x]: Port the complete active RDKit path and branched-subgraph enumeration source unit identified in Step 3, including hydrogen filtering, `fromAtoms`, duplicate/order semantics, rings, disconnected components, and source container ordering.

Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 19 [x]: Add focused enumeration tests that compare ordered bond-index paths and subgraphs against pinned RDKit for linear, branched, cyclic, fused, explicit-H, disconnected, and restricted-start fixtures.

Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 21 [x]: Run the exact focused Rust test target added in Step 19 with `op-contracts-strict`.

Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 23 [x]: Port the complete active `RDKitFPEnvGenerator` source unit identified in Step 3, including path-hash ordering, one-bond special handling, atom participation, bit paths, atom-to-bits deduplication, and atom counts.

Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 25 [x]: Add focused exact-environment tests for raw bit identifiers and every selected `AdditionalOutput` branch and interaction.

Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 27 [x]: Run the exact focused Rust test target added in Step 25 with `op-contracts-strict`.

Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 29 [x]: Port the complete active generic-generator and `RDKFingerprintMol` output assembly source unit identified in Step 3, including preconditions, Boost random-bit expansion, modulo behavior, collisions, density folding, minimum size, result resizing, and provenance copying.

Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 31 [x]: Add focused exact-bit RDKFingerprint tests for defaults, every public option, pairwise option interactions, invalid values, random-bit ordering, collisions, folding thresholds, output-size changes, and provenance after folding.

Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 33 [x]: Run the exact focused Rust RDKFingerprint source-parity test target added in Step 31 with `op-contracts-strict`.

Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 35 [x]: Replace the Rust and Python RDKFingerprint fail-closed dispatch with the source-backed implementation without adding compatibility shims for removed placeholder-only parameters.

Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 37 [x]: Add focused Python API tests for RDKFingerprint defaults, keywords, typed errors, output ordering, provenance serialization, and scalar-versus-repeated-call determinism.

Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 39 [x]: Run the exact focused Python RDKFingerprint API test target added in Step 37.

Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 41 [x]: Add a pinned RDKFingerprint oracle generator and metadata profile covering all focused branches plus a source-meaningful parameter matrix over the maintained 5,000-row corpus.

Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 43 [x]: Run the RDKFingerprint oracle generator from Step 41 and preserve checksummed expected data and generation logs.

Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 45 [x]: Add an exact-bit RDKFingerprint 5,000-row parity test that consumes every generated profile without filtering, tolerance, fallback, or sampling.

Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 47 [x]: Run the exact RDKFingerprint 5,000-row parity test added in Step 45 in release mode with deterministic parallel execution.

Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 49 [x]: Audit RDKFingerprint source-line closure, exact parity, error parity, determinism, allocation shape, asymptotic complexity, and benchmark results and write `dev/gap_reports/rdkit_topological_fingerprint_full_port_validation.md` with links to ignored machine artifacts.

Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 51 [x]: Audit RDKit's Avalon adapter source, `ava-formake` tag and checksum, every reachable source-file license, redistribution obligations, public defaults, and current COSMolKit placeholder mismatch and write `dev/gap_reports/avalon_fingerprint_source_provenance.md`.

Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 53 [x]: Materialize the checksum-verified `AvalonToolkit_2.0.5-pre.3` archive under ignored `tmp/parity-audit/sources/avalon/` according to the provenance report from Step 51.

Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 55 [x]: Audit the complete `getAvalonFP(ROMol, ...)` reachable call graph across the RDKit adapter and Avalon C engine and write `dev/gap_reports/avalon_fingerprint_source_inventory.md` with exact source ranges, flag-family ownership, prerequisite closure, a line-coverage ledger, and a public-option mapping.

Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 57 [x]: Modify the still-unsupported Rust and Python Avalon parameter types to the selected source boundary, replace invented path options with source-backed fields, encode bit flags as a typed mask, and retain structured fail-closed behavior.

Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 59 [x]: Add focused Rust tests for corrected Avalon defaults, Python-versus-C++ default profiles, typed flag validation, output-size boundaries, and structured unsupported results.

Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 61 [x]: Run the exact focused Rust Avalon contract test target added in Step 59 with `op-contracts-strict`.

Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 63 [x]: Port the complete active RDKit `molToReaccs` conversion and Avalon REACCS data/parser source unit identified in Step 55, including MolBlock serialization, locale-independent numeric parsing, atom and bond fields, query fields, and allocation/error behavior.

Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 65 [x]: Add focused full-state conversion tests comparing COSMolKit's internal Avalon molecule state with an instrumented pinned native adapter for ordinary, aromatic, isotope, charge, radical, query, explicit-H, and boundary-size fixtures.

Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 67 [x]: Run the exact focused Avalon conversion test target added in Step 65 with `op-contracts-strict`.

Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69a [x]: Port the complete Avalon `SetupNeighbourhood` source unit used by fingerprint preprocessing, including source bond-table ordering, atom/bond index pairing, the source `MAXNEIGHBOURS` boundary, and structured handling before the C fixed-array overflow boundary.

Step 69b [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69c [x]: Port the complete Avalon fingerprint hydrogen-count closure, including `GuessHCountsFromSubstitution`, `ComputeImplicitH`, `ImplicitHydrogens`, the full source valence table, query-H initialization, explicit-isotope hydrogen handling, and source one-based array semantics.

Step 69d [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69e [x]: Port the complete Avalon basis-ring closure used by fingerprint preprocessing, including `RingList`, `SortRings`, `CombineRings`, bit-set operations, bond-order-dependent spanning-tree order, and the pinned C-runtime random sequence used for equal-cardinality toggles.

Step 69e.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69e.2 [x]: Port the complete Avalon fingerprint `RingState`, `MarkRecursive`, and `SetRingSizeFlags` closure, including atom/bond ring status, recursive ring-size traversal, source bit flags, and all restoration behavior local to those routines.

Step 69f [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69g [x]: Port the complete Avalon `ProperRingPairs` helper and `PerceiveAromaticBonds` fingerprint aromaticity mode, including basis-plus-proper-ring ordering, cumulene exclusion, repeated propagation, and exact bond mutation order.

Step 69g.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69g.2 [x]: Port the complete Avalon `AtomSymbolMatch` helper and `PerceiveDYAromaticity` symbol-list lookup closure, including every source pseudo-symbol table and exact comma-token and whole-pattern fallback order.

Step 69g.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69g.4 [x]: Port the complete Avalon `PerceiveDYAromaticity` candidate-ring construction unit, including initial ring membership, candidate-subgraph reconstruction, basis reduction, proper-ring insertion, fused-pair generation, and exact linked-list order.

Step 69g.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69g.6 [x]: Port the remaining complete Avalon `PerceiveDYAromaticity` evaluation and mutation unit, including local pi-electron rules, exocyclic pull, charge exclusions, query-H mutation, repeated propagation, exact bond mutation order, and final aromatic single-bond closure.

Step 69h [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69i [x]: Port the complete active Avalon atom/bond preprocessing color and counter closure, including degree, carbon degree, unsaturation, special-neighbour, fusion-state, symbol-list, and query-mode branches.

Step 69j [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69k [x]: Port the complete shared Avalon fingerprint traversal and hashing closure, including path recursion, special-neighbour recursion, feature-pair distances, ring-size pairs, hash arithmetic, and source mutation restoration around the complete preprocessing lifecycle.

Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 71 [x]: Add focused intermediate-state tests for every preprocessing array, traversal order, hash boundary, aromaticity pass, ring annotation, and restoration path used by Avalon fingerprint generation.

Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 73 [x]: Run the exact focused Avalon preprocessing test target added in Step 71 with `op-contracts-strict`.

Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 75 [x]: Port the complete Avalon flag-family source units `USE_RING_PATTERN` through `USE_ATOM_COUNT` identified in Step 55 with verbatim line-by-line two-axis source anchors.

Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 77 [x]: Add focused exact-count and exact-bit tests for each Avalon flag from `0x000001` through `0x000010`, every pairwise interaction in that group, query mode, and both aromaticity passes.

Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 79 [x]: Run the exact focused Avalon low-flag test target added in Step 77 with `op-contracts-strict`.

Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 81 [x]: Port the complete Avalon flag-family source units `USE_AUGMENTED_ATOM` through `USE_AUGMENTED_BOND` identified in Step 55 with verbatim line-by-line two-axis source anchors.

Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 83 [x]: Add focused exact-count and exact-bit tests for each Avalon flag from `0x000020` through `0x000400`, every pairwise interaction in that group, explicit/implicit hydrogen distinctions, query mode, and both aromaticity passes.

Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 85 [x]: Run the exact focused Avalon middle-flag test target added in Step 83 with `op-contracts-strict`.

Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 87 [x]: Port the complete Avalon flag-family source units `USE_RING_SIZE_COUNTS` through `USE_FEATURE_PAIRS` identified in Step 55 with verbatim line-by-line two-axis source anchors.

Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 89 [x]: Add focused exact-count and exact-bit tests for each Avalon flag from `0x000800` through `0x004000`, every pairwise interaction in that group, fused/spiro rings, class spiders, feature-pair distances, and source boundary sizes.

Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 91 [x]: Run the exact focused Avalon high-flag test target added in Step 89 with `op-contracts-strict`.

Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 93 [x]: Port the complete Avalon non-SSS flag-family source units `USE_SCAFFOLD_IDS`, `USE_SCAFFOLD_COLORS`, `USE_SCAFFOLD_LINKS`, and `USE_SHORTCUT_LABELS` identified in Step 55 with verbatim line-by-line two-axis source anchors.

Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 95 [x]: Add focused exact-count and exact-bit tests for every Avalon non-SSS flag, their combined `0xF00000` profile, shortcut labels, scaffold topology, exclusion behavior, and non-query gating.

Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 97 [x]: Run the exact focused Avalon non-SSS test target added in Step 95 with `op-contracts-strict`.

Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 99 [x]: Port the complete active Avalon `SetFingerprintCountsWithFocus`, `SetFingerprintBits`, RDKit adapter, and explicit-bit-vector assembly source unit identified in Step 55, including `isQuery`, the second Daylight-aromaticity accumulation pass, `resetVect`, `bitFlags`, `nBits / 8`, four-byte internal rounding, signed-char promotion, partial-byte behavior, and cleanup paths.

Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 101 [x]: Add focused exact-bit Avalon wrapper tests for every public option, default profile, all-flags profile, zero and invalid sizes, non-multiples of eight, non-multiples of thirty-two, query mode, reset/accumulate internals, repeated calls, and allocation/error paths.

Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 103 [x]: Run the exact focused Rust Avalon source-parity test target added in Step 101 with `op-contracts-strict`.

Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 105 [x]: Replace the Rust and Python Avalon fail-closed dispatch with the source-backed implementation without compatibility shims for removed placeholder-only path parameters.

Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 107 [x]: Add focused Python Avalon API tests for defaults, typed flags, keyword interactions, query mode, output-size boundaries, typed errors, repeated-call determinism, and non-mutation of the source molecule.

Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 109 [x]: Run the exact focused Python Avalon API test target added in Step 107.

Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 111 [x]: Add a pinned Avalon oracle generator and metadata profile covering every flag family, query mode, size boundary, aromaticity pass, focused regression case, and a source-meaningful parameter matrix over the maintained 5,000-row corpus.

Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 113 [x]: Run the Avalon oracle generator from Step 111 and preserve checksummed expected data and generation logs.

Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 115 [x]: Add an exact-bit Avalon 5,000-row parity test that consumes every generated profile without filtering, tolerance, fallback, or sampling.

Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 117 [x]: Run the exact Avalon 5,000-row parity test added in Step 115 in release mode with deterministic parallel execution.

Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 119 [x]: Audit Avalon source-line closure, exact parity, error parity, determinism, allocation shape, asymptotic complexity, native-adapter performance, and license-notice completeness and write `dev/gap_reports/avalon_fingerprint_full_port_validation.md` with links to ignored machine artifacts.

Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 121 [x]: Update family-specific feature metadata, support inventory, parity scope, Rust docs, Python docs, generated stubs, README, examples, and release-facing unsupported text only for families whose full-port validation report has zero unexplained mismatches.

Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 123 [x]: Add combined Rust and Python tests for calling RDKFingerprint, Avalon, Morgan, and MACCS sequentially and concurrently on the same molecules without state mutation, cache leakage, option cross-talk, ordering drift, or thread nondeterminism.

Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 125 [x]: Run the exact combined fingerprint interaction test targets added in Step 123.

Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 127 [x]: Run formatting, strict workspace checks, full release workspace tests, Python tests, generated-stub verification, rustdoc, Sphinx warnings-as-errors, fingerprint examples, feature-combination checks, exact oracle reproducibility, and deterministic performance comparisons as the final validation matrix.

Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 129 [x]: Write `dev/gap_reports/rdkit_topological_avalon_fingerprint_completion_audit.md` proving every scope item, source line, option branch, corpus row, public document, test gate, performance requirement, and support-status change from this plan is complete or identifying the exact blocking evidence.
