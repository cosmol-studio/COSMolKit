# RDKit MolAlign and RMSD Full Source-Port Plan

## Execution Contract

- This plan is for sequential continuous execution.
- Execute unchecked steps in order and continue until every step is completed, blocked, or interrupted.
- Mark a completed step by changing only its `[ ]` to `[x]`.
- Every real task step must be immediately preceded by a fresh reading step for `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
- Reading steps must be executed, not treated as already read.
- `Audit`, `Port`, `Implement`, `Modify`, `Update`, and `Fix` steps must produce a concrete artifact.
- A test-writing step must be followed immediately by its most specific test command after the required reading step.
- Final whole-plan validation is required after all implementation steps.
- Do not add heuristic chemistry, silent fallbacks, placeholder branches, or a second implementation of an existing capability.
- During porting, do not execute git operations unless explicitly requested by the human author.

## Scope And Boundary

- The pinned source is RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`, under `third_party/rdkit/`.
- The in-scope source boundary is `Code/Numerics/Alignment/AlignPoints.{h,cpp}`, `Code/GraphMol/MolAlign/AlignMolecules.{h,cpp}`, the relevant `MolAlign` Python wrapper overloads, `MolTransforms::transformConformer`, and the source-reachable substructure and query helpers used by those functions.
- O3A, Crippen/O3A scoring orchestration, `O3AAlignMolecules.{h,cpp}`, and its FFI/backend dependencies are explicitly out of scope and must be recorded as a separate future capability.
- RMSD calculation APIs are read-only by default; no function named as a calculation, measurement, transform query, or `rmsd` may mutate a molecule.
- Any coordinate mutation must use an explicit project-native mutating name with a trailing `_`, or a value-style method returning a new molecule/result; RDKit's implicit `getBestRMS` probe mutation must not be reproduced under an ambiguous calculation name.

## Backend Policy

- Production `cosmolkit-core` has exactly one private source-equivalent `Numerics::Alignment` implementation for quaternion construction, Jacobi diagonalization, reflection, weighted/unweighted sums, SSR, translation, and transform application.
- Existing `coordinates.rs` and `distgeom.rs` alignment implementations must be consolidated into that core without changing their validated floating-point operation order or behavior.
- The pinned `AlignPoints.cpp` path has no external mathematical FFI dependency, so this task must not manufacture an FFI production dependency merely to satisfy a backend shape.
- Introduce FFI only if exact parity requires a numerical behavior that Rust cannot reproduce directly; source availability or convenience alone is not sufficient justification.
- The audited `AlignPoints.cpp` boundary is fully expressible with source-ordered Rust floating-point operations, so ordinary MolAlign has one pure-Rust production implementation and pinned-RDKit generated oracle data, with no FFI backend feature.
- Backend selection, numerical tolerance, source version, and compared observable fields must be recorded in fixtures, manifests, and parity reports.

## Normative Call Graph

```text
Rust/Python MolAlign API
  -> typed parameter validation and conformer-ID resolution
  -> one private Numerics::Alignment core
  -> explicit atom-map validation or canonical RDKit substructure mapping
  -> optional terminal conjugated-group symmetry query on a match-only clone
  -> weighted/unweighted transform or coordinate-frame RMS calculation
  -> read-only result, or explicitly named coordinate mutation operation
```

## Function-By-Function Source Ledger

| Source unit | Required port/reuse action |
|---|---|
| `RDNumeric::Alignments::AlignPoints` | Consolidate into the sole private Rust alignment core; preserve weights, positive-weight preconditions, reflection, max iterations, tolerance, SSR clamp, translation order, and empty/mismatched input errors. |
| `RDNumeric::Alignments::jacobi` | Reuse the existing source-marked Jacobi implementation after operation-order audit; remove duplicate copies only after focused parity proves unchanged results. |
| `_weightedSumOfPoints`, `_weightedSumOfLenSq`, `_sumOfWeights`, `_computeCovarianceMat`, `_covertCovMatToQuad` | Consolidate typed fixed-array helpers in the same core, including weighted and unweighted paths and 32/64-bit-independent arithmetic semantics. |
| `RDGeom::Transform3D` construction/application | Reuse one fixed-array transform representation and one point-application helper for alignment, depiction, distgeom, and MolAlign. |
| `alignConfsOnAtomMap` | Port as a private map-to-point adapter with exact probe/reference orientation and RMS normalization. |
| `calcMSDInternal` | Port coordinate-frame weighted RMS calculation without alignment or mutation. |
| `getAllMatchesPrbRef` | Reuse the canonical substructure matcher with `uniquify=false`, recursive queries, source-specific query-query behavior, exact max-match handling, and optional match-only terminal-group symmetrization. |
| `symmetrizeTerminalAtoms` | Promote the existing distgeom helper to the sole shared match-query helper; preserve clone-only topology/query edits and cache lifecycle. |
| `getBestRMSInternal` | Port exact candidate iteration, strict `<` tie selection, optional transform/map output, thread partition/merge ordering, and structured no-match/weight errors. |
| `getAlignmentTransform` | Add explicit-map and source-backed automatic-map paths with RDKit defaults, conformer IDs, transform result, and no mutation. |
| `alignMol` | Expose only as an explicitly named value-style or trailing-underscore coordinate mutation operation; delegate transform calculation to the shared core and preserve operation-contract metadata. |
| `getBestAlignmentTransform` | Return best RMS, transform, and selected map without changing either molecule. |
| `getBestRMS` | Provide a project-native read-only best-RMS API; if an in-place aligned-probe variant is exposed, give it an unmistakable trailing-underscore mutation name and separate result type. |
| `getAllConformerBestRMS` | Port triangular pair ordering, real conformer-ID lookup, mapping reuse, symmetry options, and deterministic threaded reduction. |
| `CalcRMS` | Port no-alignment coordinate-frame RMS, map enumeration, symmetry option, weights, and no mutation. |
| `alignMolConformers` | Port selected/all conformer behavior, atom selection, weights, reflection, max iterations, first-reference semantics, RMS list ordering, and explicit coordinate mutation naming. |
| Python wrapper overloads | Add only project-native names and conversions; document any intentional divergence from RDKit's mutating `getBestRMS` surface. |

## Public API And Contract Design

- Define typed alignment parameters, atom-map pairs, transform/result types, conformer selectors, and structured alignment errors in the owning core module.
- Resolve conformers by RDKit ID semantics (`-1` means first; nonnegative values match stored IDs), then use one resolver everywhere in MolTransforms, MolAlign, and tests.
- Keep query/topology reads narrow and immutable; automatic mapping may use the existing matcher but must not expose mutable molecule storage.
- Register every coordinate-mutating operation through `molecule_ops!` with coordinate write access, explicit `may_mutate`, drawing/cache effects, preservation of topology/CIP state, and invariant/parity metadata.
- Do not register read-only RMSD or transform queries as mutation operations.
- Keep one core implementation for scalar and batch APIs; wrappers must delegate rather than reimplement candidate enumeration or numerical math.

## Required Test And Oracle Matrix

- Alignment-core tests cover empty, one-point, collinear, coincident, weighted, zero/negative weight, mismatched length, reflection, max-iteration, near-zero negative SSR, transform translation, and exact operation ordering.
- Consolidation tests prove coordinates, distgeom pruning, and MolAlign all use the same core and preserve their existing outputs.
- Conformer tests cover first/real ID selection, sparse/non-sequential IDs, invalid IDs, multiple conformers, selected subsets, and triangular output order.
- Atom-map tests cover explicit maps, automatic probe-to-reference orientation, duplicated/reordered maps, no-match errors, max-match limits, hydrogen filtering, and exact tie selection.
- Symmetry tests cover nitro/carboxylate-like terminal groups, no-match molecules, query-only clone behavior, aromaticity/charge state, and cache preservation/invalidation.
- Public API tests cover Rust and Python read-only immutability, explicit mutating names, result transforms/maps/RMS values, errors, weights, reflection, max iterations, and concurrency.
- Pinned-RDKit oracle tests compare the complete transform matrix, RMS/SSR, selected map, candidate ordering, and errors through generated expected data without linking FFI into the Rust test or production dependency graph.
- Pure-Rust parity tests compare against pinned RDKit with documented numerical tolerances and fail closed on unsupported source branches. There is no backend selection while Rust can reproduce the source behavior.
- Maintained corpus gates cover the curated alignment corpus, the 5,000-row corpus, and mutually parseable ChEMBL 37 records where 3D conformers and maps are available; report attempted, parsed, compared, unsupported, errored, and mismatching counts without filtering failing rows.
- Regression fixtures lock every discovered discrepancy before a fix is accepted.

## Completion Criteria

- All in-scope ledger rows have source-marked code or an explicit exact-reuse proof; no duplicate alignment implementation remains.
- `coordinates`, `distgeom`, and MolAlign delegate to one private `Numerics::Alignment` core.
- Pure-Rust production behavior is source-backed and no FFI is introduced unless a future measured parity blocker proves Rust cannot reproduce a required numerical operation.
- Read-only RMSD and transform APIs never mutate molecule state; every mutating path is explicit in its name and operation contract.
- Conformer ID, mapping direction, weights, symmetry, threading, errors, and observable result state match the pinned source boundary.
- Focused, 5,000-row, and applicable ChEMBL 37 parity gates pass with no unexplained mismatch.
- Rust/Python exports, docs, examples, validation scope, inventory, and changelog describe ordinary MolAlign only and keep O3A as a separate future capability.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit the pinned AlignPoints and MolAlign call graph and write `dev/gap_reports/rdkit_molalign_source_inventory.md` with exact source ranges, observable behavior, dependencies, and O3A boundary.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit existing coordinate, distgeom, transform, conformer-ID, matcher, symmetry, operation, Python, and test infrastructure and append concrete reuse and duplication findings to `dev/gap_reports/rdkit_molalign_source_inventory.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Define the typed MolAlign parameter, result, transform, conformer-selector, atom-map, error, and numerical-tolerance contracts in the source inventory and owning API design artifact.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Add focused contract tests for the typed MolAlign API boundary, read-only immutability, conformer-ID semantics, and explicit mutation naming.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Run the focused API-boundary tests added in Step 9 with `cargo test -p cosmolkit-core --release --features op-contracts-strict` and the corresponding Python filter when available.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Consolidate the existing AlignPoints math into the sole private `Numerics::Alignment` core with complete source markers and preserve validated callers.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Add focused alignment-core parity and error tests for weighted/unweighted, reflection, Jacobi, transforms, SSR, and iteration behavior.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Run the focused alignment-core tests added in Step 15 with `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Consolidate conformer-ID resolution and the terminal-group symmetrization helper into shared private source-backed adapters used by all callers.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Add focused conformer-ID and terminal-symmetry tests covering source defaults, sparse IDs, clone-only edits, and cache behavior.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Run the focused conformer and terminal-symmetry tests added in Step 21 with `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Port explicit-map `getAlignmentTransform`, coordinate-frame `CalcRMS`, and their complete structured validation behavior through the shared alignment core without mutation.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Add focused explicit-map transform and CalcRMS parity tests for values, errors, weights, reflection, IDs, and immutability.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run the focused explicit-map and CalcRMS tests added in Step 27 with `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Port automatic mapping, candidate enumeration, BestAlignmentTransform, BestRMS, and exact tie/thread/max-match behavior through the shared matcher and alignment core.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Add focused automatic-map and BestRMS parity tests for symmetry, hydrogens, duplicate maps, tie order, no matches, and complete observable results.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Run the focused automatic-map and BestRMS tests added in Step 33 with `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Port all-conformer BestRMS and align-conformers behavior with triangular ordering, selected IDs/atoms, weights, reflection, and explicit mutation semantics.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Add focused all-conformer and align-conformers tests for ordering, selection, coordinate changes, operation effects, errors, and read-only separation.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Run the focused all-conformer and align-conformers tests added in Step 39 with `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Regenerate the remaining plan steps to include typed operation output, immediate focused validation, and the current root `VALIDATION.md` path.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Extend `molecule_ops!` with an optional typed operation output and make conformer alignment return its source-produced ordered RMS report from the same coordinate transaction.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Add focused macro and public API tests for typed value and in-place output, successful commit, structured errors, source preservation, and the absence of duplicate alignment calculation.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run the focused typed-operation and MolAlign tests with `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Add pinned-RDKit focused expected-data generation, complete-state oracle tests, source-version checks, and dependency guards enforcing the sole pure-Rust runtime path.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Run the focused pinned-RDKit MolAlign oracle and pure-Rust dependency-boundary tests.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Add project-native Python MolAlign parameter, result, read-only, value-mutation, and trailing-underscore APIs with generated stubs, examples, and binding tests.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run the focused Python MolAlign tests, regenerate stubs, and run the relevant type checks.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Add maintained curated and 5,000-row RDKit MolAlign fixtures with complete transform, RMSD, map, conformer, error, and mutation-state schemas.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run the maintained curated and 5,000-row MolAlign parity gates with release strict contracts.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Add the reproducible ChEMBL 37 MolAlign audit profile and repository-owned preparation and run manifests for supported 3D inputs.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Run the complete applicable ChEMBL 37 MolAlign profile and record exact attempted, parsed, compared, errored, unsupported, and mismatch counts.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Audit source markers, duplicate implementations, public exports, operation matrices, the pure-Rust numerical boundary, and the separate unsupported O3A boundary and write the final gap report.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Update `VALIDATION.md`, `dev/porting_inventory.md`, public Rust and Python documentation, examples, API design, and the unreleased changelog with the validated ordinary MolAlign boundary and separate O3A scope.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Run `cargo fmt --all`, strict core check and release tests, strict workspace release tests, Python tests, generated-stub checks, type checks, and the Sphinx build as the final whole-plan gate.
