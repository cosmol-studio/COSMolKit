# RDKit Pattern Fingerprint Full Source-Port Plan

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
- Never treat a required reading step as "already read".
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
- Do not use "smallest subpart", skeleton code, dispatch-only code, placeholder code, TODO-only code, or partial porting as completion.
- If a step is too large to complete fully, regenerate the plan by splitting that step into smaller full-port steps.
- Each split step must still represent a complete source-backed behavior, not a placeholder.
- A `Port` step is complete only when the selected RDKit behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Source Pin And Boundary

- Oracle and source tree: RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`, vendored under `third_party/rdkit/`.
- Source root: `Code/GraphMol/Fingerprints/PatternFingerprints.cpp::PatternFingerprintMol(const ROMol&)` and its reachable SMARTS, query, ring, substructure, hash, bit-vector, and Python-wrapper functions.
- The complete ordinary-molecule boundary includes all 13 active built-in SMARTS patterns, `fpSize`, ordinary and tautomeric modes, and ordinary/query-bearing molecules.
- In the pinned ordinary overload, `atomCounts` and `setOnlyBits` are size-checked but otherwise unused; exact porting preserves this observable source behavior instead of implementing the header comment's intended behavior.
- `PatternFingerprintMol(const MolBundle&)` is excluded until COSMolKit deliberately defines a public bundle model; its member-fingerprint intersection rule remains documented source evidence, not a disguised slice overload.
- Public APIs remain project-native and method-oriented and return the existing `Fingerprint` representation.

## Normative Call Graph

```text
Rust/Python Pattern API
  -> PatternFingerprintParams validation
  -> pattern_fingerprint_core
     -> compiled_pattern_fingerprint_queries
        -> existing mol_from_smarts for the 13 pinned patterns
     -> existing fast_find_rings when ring information is below fast-or-better
     -> source-backed atom-query complexity plus Pattern AtomNull override
     -> isPatternComplexQuery
     -> isTautomerBondQuery
     -> existing substructure matcher
        -> uniquify=false
        -> maxMatches=100000000
        -> remaining pinned defaults
     -> updatePatternFingerprint
        -> existing 32-bit gboost hash_combine
        -> existing Fingerprint bit insertion
  -> scalar and ordered batch facades delegate to the same core
```

## Function-By-Function Source Ledger

| Source function or unit | Pinned source | Required action |
|---|---|---|
| `isComplexQuery(const Bond*)` | `Code/GraphMol/QueryOps.cpp:728` | Consolidate exact no-query, negation, description, `BondOr`, and Single-or-Aromatic handling into the sole crate-private classifier. |
| `_complexQueryHelper` | `Code/GraphMol/QueryOps.cpp:769` | Consolidate recursive atom-query traversal with sibling-carried `hasAtNum`. |
| `isComplexQuery(const Atom*)` | `Code/GraphMol/QueryOps.cpp:892` | Consolidate AtomNull, atomic-number/type, Or/Xor, And, negation, and missing-atomic-number behavior. |
| `MolOps::fastFindRings` | `Code/GraphMol/FindRings.cpp:1131` | Reuse the existing complete port; do not add Pattern-local ring perception. |
| `SmartsToMol` | `Code/GraphMol/SmilesParse/SmilesParse.h:189` | Reuse the completed SMARTS port as the only built-in pattern compiler. |
| `SubstructMatch` vector overload | `Code/GraphMol/Substruct/SubstructMatch.cpp:481` | Reuse the completed matcher with exactly the two Pattern parameter overrides. |
| `gboost::hash_combine` | `Code/RDGeneral/hash/hash.hpp:211` | Reuse the existing sole wrapping 32-bit implementation. |
| `ExplicitBitVect::setBit` | `Code/DataStructs/ExplicitBitVect.cpp:84` | Reuse existing `Fingerprint` insertion and prove collisions and repeated writes. |
| `ss_matcher::ss_matcher()` | `PatternFingerprints.cpp:41` | Omit if the Rust cache type needs no meaningless default constructor. |
| `ss_matcher::ss_matcher(const std::string&)` | `PatternFingerprints.cpp:42` | Compile repository-owned constants through the existing SMARTS parser and treat failure as an internal invariant. |
| `ss_matcher::getMatcher` | `PatternFingerprints.cpp:49` | Return immutable cached query references without per-molecule clones. |
| `pqs` and `pattern_flyweight` | `PatternFingerprints.cpp:57-80` | Port the exact 13 active strings and lazy thread-safe no-eviction caching semantics. |
| `detail::getAtomNumbers` | `PatternFingerprints.cpp:83` | Record as unreachable from both pinned Pattern overloads and do not add dead code. |
| `isPatternComplexQuery` | `PatternFingerprints.cpp:140` | Port exact no-query, negation, and description-not-`BondOrder` behavior. |
| `isTautomerBondQuery` | `PatternFingerprints.cpp:154` | Port exact recognition of the two tautomer bond-query descriptions. |
| `updatePatternFingerprint` | `PatternFingerprints.cpp:161` | Port preconditions, cache iteration, ring preparation, masks, matching, occurrence hashing, atom/bond hashing, query suppression, tautomer hashing, collisions, and inert arguments. |
| `PatternFingerprintMol(const ROMol&)` | `PatternFingerprints.cpp:339` | Implement the ordinary-molecule facade as the sole Pattern production entry. |
| `PatternFingerprintMol(const MolBundle&)` | `PatternFingerprints.cpp:355` | Preserve an explicit out-of-scope proof pending a public bundle model. |
| `wrapPatternFingerprint` | `Code/GraphMol/Wrap/MolOps.cpp:666` | Use as defaults, width, conversion, mutable-list, and exception oracle while retaining project-native names. |
| Version constant | `Fingerprints.h:149` | Preserve source version `1.0.0` as read-only metadata. |

## Required Test Matrix

- Query classification covers absent, null, atomic-number/type, aromatic/aliphatic, And/Or/Xor, recursion, negation, Single-or-Aromatic, and Single-or-Double-or-Aromatic query trees.
- Molecules cover empty and isolated atoms, disconnected components, chains, branches, rings of varied size, fused/bridged/spiro systems, aromatic/kekulized pairs, hydrogens, charges, isotopes, unusual bonds, and query-bearing graphs.
- Every built-in pattern is exercised with zero, one, repeated, overlapping, and symmetry-related matches using exact non-unique match multiplicity and order.
- Intermediate behavior covers occurrence hash evolution, atom-map order, atomic-number hashes, bond-order hashes, aromatic normalization, query suppression, tautomer-query suppression, `u32::MAX` tautomer hashing, modulo, collisions, and single-atom patterns.
- Options cover zero and boundary/non-power-of-two widths, tautomer mode, correctly/incorrectly sized count/mask inputs, pre-populated counts, and all-zero/all-one/sparse masks, including exact proof that valid ordinary-overload count/mask values are inert.
- Screening tests require exact vectors and all-probe-bits containment for every pinned upstream and project query/target pair where the reference matcher reports a match.
- API tests cover Rust/Python scalar and ordered batch calls, errors, immutability, repeated calls, shared values, concurrency, and interleaving with every implemented fingerprint family.
- Corpus gates cover every active focused row, `smiles_small`, `smiles_5000`, and every mutually parseable ChEMBL 37 row for every active profile without sampling, filtering, allowlists, tolerance, or fallback.
- Large-corpus comparison records exact width and complete ordered on-bit set plus attempted, parsed, compared, errored, and mismatching row counts.

## Completion Criteria

- Every ledger row has a concrete source-marked code/test artifact or exact unreachable/out-of-scope proof.
- There is one query model, SMARTS compiler, matcher, ring engine, hash implementation, vector representation, Pattern table/cache, and Pattern core.
- Pattern compilation is lazy, thread-safe, shared, and performs no per-molecule query clone or reparse.
- Ordinary, tautomeric, and query-bearing behavior including pinned inert arguments matches exactly.
- Focused, small, 5,000-row, and complete ChEMBL 37 profiles have zero unexplained mismatch with reproducible identities and exact counts.
- Benchmarks show no duplicate compilation, avoidable molecule cloning, worse matcher complexity, or unexplained hot-path regression.
- Rust/Python exports, stubs, docs, examples, metadata, inventory, parity scope, and changelog describe only the validated boundary.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit the complete pinned ordinary-molecule Pattern call graph and write `dev/gap_reports/rdkit_pattern_fingerprint_source_inventory.md` with exact ranges, line coverage, reachable/dead helpers, overload boundary, defaults, errors, upstream tests, and dependency ownership.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit current SMARTS, substructure, ring, query, hash, vector, fingerprint, batch, Python, oracle, and ChEMBL facilities and append a concrete reuse and duplication gap report to `dev/gap_reports/rdkit_pattern_fingerprint_source_inventory.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Consolidate the complete reachable `_complexQueryHelper` and atom/bond `isComplexQuery` units into the sole narrow crate-private typed helpers with line-level markers and exhaustive local tests.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run the exact focused Rust query-complexity test filter added in Step 7 with `op-contracts-strict`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port the exact Pattern built-in table and matcher flyweight as one lazy thread-safe immutable compiled-query set with constant-order, compile-once, concurrent-access, and no-clone tests.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Run the exact focused Rust Pattern compiled-query test target added in Step 11 with `op-contracts-strict`.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Port `isPatternComplexQuery` and `isTautomerBondQuery` as typed query-tree classifiers with exact description and negation tests.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Run the exact focused Rust Pattern query-classifier test filter added in Step 15 with `op-contracts-strict`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port the complete `updatePatternFingerprint` function for ordinary and query-bearing molecules with line-level markers and intermediate tests for every pattern, match-count hash, atom hash, bond hash, query branch, tautomer branch, collision, and pinned inert argument.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run the exact focused Rust Pattern core test target added in Step 19 with `op-contracts-strict`.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Port the ordinary-molecule Pattern facade into one project-native Rust parameter/result API with exact defaults, widths, validation, metadata, and focused public tests.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Run the exact focused Rust Pattern public API test target added in Step 23 with `op-contracts-strict`.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Add project-native Python Pattern scalar and ordered batch methods delegating to the sole Rust core with type, default, error, immutability, ordering, and concurrency tests.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run the exact focused Python Pattern API test target added in Step 27.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Add pinned Pattern profiles, oracle generator, expected-data schema, version assertion, identity manifest, and focused/small/5,000-row registration under the existing testdata workflow.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run Pattern golden preparation for every focused, small, and 5,000-row profile and verify all manifests and checksums.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Add exact Pattern focused, complete-small, and complete-5,000-row Rust parity tests for every active ordinary, tautomeric, query, size, count-state, and mask-validation profile without row filtering.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run the exact Pattern focused, complete-small, and complete-5,000-row parity targets added in Step 35 in release mode with `op-contracts-strict`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Add the complete pinned upstream Pattern screening regressions and project query/target branch matrix as exact-vector and all-probe-bits tests.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Run the exact focused Rust Pattern screening test target added in Step 39 with `op-contracts-strict`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Extend the repository-owned ChEMBL 37 workflow with complete deterministic Pattern profiles, exact schemas, sharding, identity, resumability, and zero-mismatch acceptance tests.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run the exact focused workflow tests for the Pattern ChEMBL additions from Step 43.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Run every active Pattern profile over every mutually parseable ChEMBL 37 row with deterministic parallel execution and preserve the complete external run identity, counts, checksums, and mismatch artifacts.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Add representative Pattern benchmarks against pinned RDKit for path-rich, ring-rich, query-bearing, disconnected, and large molecules with allocation and compiled-cache reuse measurements.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run the Pattern benchmark suite from Step 49 and preserve its environment, raw results, and complexity review outside tracked machine artifacts.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Add Pattern composition tests for repeated, interleaved, and concurrent calls with every implemented fingerprint family on shared molecules without mutation, cache leakage, ordering drift, or option cross-talk.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Run the exact Rust and Python Pattern composition targets added in Step 53 with `op-contracts-strict` on Rust.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Add static architecture tests rejecting duplicate Pattern tables, compilation loops, query classifiers, match loops, hash implementations, vectors, and scalar/batch cores.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Run the exact static Pattern architecture test target added in Step 57 with `op-contracts-strict`.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run `cargo fmt --all`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Regenerate `python/cosmolkit.pyi` with the project stub generator.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Add focused stub-surface tests for Pattern scalar and batch signatures.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Run the exact Pattern stub-surface tests added in Step 69.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run `.venv/bin/basedpyright python/tests python/examples`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run the focused Pattern Python test suite.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Audit Pattern source-line, branch, option, output, error, allocation, complexity, API, and corpus closure and write `dev/gap_reports/rdkit_pattern_fingerprint_full_port_validation.md` with exact evidence and no incomplete success claim.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Update Pattern Rust/Python docs and examples with exact options, upstream experimental metadata, query behavior, inert arguments, errors, composition, and `MolBundle` exclusion.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run the repository documentation build, link, and doctest validation entrypoint for Pattern pages and examples.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Update Pattern feature metadata, inventory, parity scope, README status, and current Unreleased changelog from the completed validation report without modifying historical plans.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Run repository documentation and status consistency checks for the Pattern claims.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [ ]: Re-run every focused and complete Pattern parity target after documentation and generated-artifact changes.
Step 88 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [ ]: Run the final strict workspace release suite after all Pattern artifacts are complete.
