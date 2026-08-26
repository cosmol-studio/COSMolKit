# RDKit Layered Fingerprint Full Source-Port Plan

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
- Source root: `Code/GraphMol/Fingerprints/Fingerprints.cpp::LayeredFingerprintMol()` and its reachable QueryOps, ring, path/subgraph, hash, bit-vector, and Python-wrapper functions.
- The complete boundary includes all six active layers, arbitrary `layerFlags`, `minPath`, `maxPath`, `fpSize`, pre-populated `atomCounts`, `setOnlyBits`, linear/branched enumeration, and ordered `fromAtoms` behavior.
- The source declares ten layer slots but defines six active encoders; inactive high flags remain accepted and produce no additional layer instead of becoming invented features.
- Public APIs remain project-native and method-oriented and return the existing `Fingerprint`; atom counts use a typed result without exposing mutable internal storage.

## Normative Call Graph

```text
Rust/Python Layered API
  -> LayeredFingerprintParams validation
  -> layered_fingerprint_core
     -> existing find_sssr when ring information is uninitialized
     -> source-backed query complexity and query-aware aromaticity
     -> existing linear path or branched subgraph enumeration
        -> source-exact ordered fromAtoms aggregation
     -> six source layer encoders
     -> sorted per-layer components
     -> source-backed 32-bit gboost hash_range
     -> existing Fingerprint bit insertion
     -> setOnlyBits and once-per-accepted-path atomCounts
  -> scalar and ordered batch facades delegate to the same core
```

## Function-By-Function Source Ledger

| Source function or unit | Pinned source | Required action |
|---|---|---|
| `queryIsBondInRing` | `Code/GraphMol/QueryOps.h:303` | Expose one typed crate-private wrapper over `RingInfo::num_bond_rings`. |
| `queryBondMinRingSize` | `Code/GraphMol/QueryOps.h:309` | Expose one typed crate-private wrapper over `RingInfo::min_bond_ring_size`. |
| `isComplexQuery(const Bond*)` | `Code/GraphMol/QueryOps.cpp:728` | Consolidate exact no-query, negation, description, `BondOr`, and Single-or-Aromatic handling into the sole crate-private classifier. |
| `_complexQueryHelper` | `Code/GraphMol/QueryOps.cpp:769` | Consolidate recursive atom-query traversal with sibling-carried `hasAtNum`. |
| `isComplexQuery(const Atom*)` | `Code/GraphMol/QueryOps.cpp:892` | Consolidate AtomNull, atomic-number/type, Or/Xor, And, negation, and missing-atomic-number behavior. |
| `isAtomAromatic` | `Code/GraphMol/QueryOps.cpp:920` | Port the complete query-aware aromatic classification including compound query ordering. |
| `MolOps::findSSSR` | `Code/GraphMol/FindRings.cpp:769` | Reuse the existing complete port and preserve its initialized-state boundary. |
| `findAllSubgraphsOfLengthsMtoN` | `Code/GraphMol/Subgraphs/Subgraphs.cpp:347` | Reuse/consolidate the existing RDKFingerprint branched enumerator with exact rooted/unrooted ordering. |
| `findAllPathsOfLengthsMtoN` | `Code/GraphMol/Subgraphs/Subgraphs.cpp:443` | Reuse/consolidate the existing RDKFingerprint linear enumerator with exact Layered argument combinations. |
| `gboost::hash_combine` | `Code/RDGeneral/hash/hash.hpp:211` | Reuse the existing sole wrapping 32-bit implementation. |
| `gboost::hash_range` | `Code/RDGeneral/hash/hash.hpp:222` | Consolidate one ordered 32-bit implementation; pinned `hash_result_t` is `uint32_t` despite assignment to `unsigned long`. |
| `ExplicitBitVect::setBit` | `Code/DataStructs/ExplicitBitVect.cpp:84` | Reuse existing `Fingerprint` insertion and prove collisions/repeated writes. |
| Preconditions/ring preparation | `Fingerprints.cpp:252-267` | Port exact range and width checks and `findSSSR`, not fast-ring, preparation. |
| Query and atom caches | `Fingerprints.cpp:269-301` | Port bond ordering, three-bit bond/endpoint complexity mask, query-aware aromaticity, and atomic numbers. |
| Path selection | `Fingerprints.cpp:303-328` | Port null/present roots, linear/branched calls, duplicates, and prepend aggregation order through the shared enumerator. |
| Layer `0x01` | `Fingerprints.cpp:358-401` | Port pure-topology degree, bond-neighbor, endpoint ordering, modulo, and packing. |
| Layer `0x02` | `Fingerprints.cpp:402-424` | Port query suppression and aromatic/single normalization with exact packing. |
| Layer `0x04` | `Fingerprints.cpp:425-447` | Port endpoint-query suppression and atomic-number/degree canonicalization. |
| Layer `0x08` | `Fingerprints.cpp:448-453` | Port endpoint-query suppression and omission for non-ring bonds. |
| Layer `0x10` | `Fingerprints.cpp:454-458` | Port minimum-ring-size modulo including zero for non-ring bonds. |
| Layer `0x20` | `Fingerprints.cpp:459-473` | Port query-aware aromatic endpoint ordering and sparse packing. |
| Layer projection | `Fingerprints.cpp:475-514` | Port empty omission, sorting, distinct-atom/layer suffixes, exact hash, modulo, mask, and once-per-path count update. |
| `wrapLayeredFingerprint` | `Code/GraphMol/Wrap/MolOps.cpp:633` | Use as defaults, width, root conversion, mutable counts, and exception oracle while retaining project-native names. |
| Constants | `Fingerprints.h:112-114` | Preserve ten slots, version `0.7.0`, and substructure layer mask `0x07` as typed read-only metadata. |

## Required Test Matrix

- Preparation covers absent/initialized ring state, exact SSSR results, atom/bond index order, query masks `0..7`, atomic numbers, and every query-aware aromaticity branch.
- Enumeration covers empty/single-atom graphs, path boundaries, linear/branched/cyclic/disconnected graphs, absent/empty/duplicated/reordered roots, invalid indices, and exact aggregation order.
- Encoders cover each individual layer, cumulative prefixes, `0x07`, `0x3f`, `0xffffffff`, zero/high-only flags, degree/neighbor/ring modulo, atom-number modulo, endpoint permutations, unusual bonds, fused rings, and each query-suppression mask.
- Projection covers sort order, distinct atom count, one-based layer ID, 32-bit wrapping hash, non-power-of-two modulo, cross-layer collisions, repeated paths, all-zero/all-one/sparse masks, seeded counts, and one increment per path with any accepted layer.
- Molecules cover empty, disconnected, chains, branches, varied rings, fused/bridged/spiro systems, aromatic/kekulized pairs, explicit hydrogen, charge, isotope, unusual bond, and query-bearing graphs.
- API tests cover Rust/Python scalar and ordered batch calls, errors, counts, masks, roots, immutability, repeats, concurrency, and interleaving with every implemented fingerprint family.
- Corpus gates cover every focused row, `smiles_small`, `smiles_5000`, and every mutually parseable ChEMBL 37 row for every active profile without sampling, filtering, allowlists, tolerance, or fallback.
- Large-corpus comparison records exact width, complete ordered on-bit set, complete count vector, and attempted/parsed/compared/errored/mismatching row counts.

## Completion Criteria

- Every ledger row has a concrete source-marked code/test artifact or exact reuse proof.
- There is one query model, ring engine, linear/branched enumerator, hash-combine, hash-range, vector representation, Layered encoder set, and Layered core.
- All layers, flags, path modes, roots, masks, counts, query suppression, ring preparation, and hash-width behavior match exactly.
- Focused, small, 5,000-row, and complete ChEMBL 37 profiles have zero unexplained mismatch with reproducible identities and exact counts.
- Benchmarks show no duplicate path enumeration, avoidable molecule cloning, worse asymptotics, or unexplained hot-path regression.
- Rust/Python exports, stubs, docs, examples, metadata, inventory, parity scope, and changelog describe only the validated boundary.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit the complete pinned Layered call graph and write `dev/gap_reports/rdkit_layered_fingerprint_source_inventory.md` with exact ranges, line coverage, dependency/reuse boundaries, defaults, errors, upstream tests, and dependency ownership.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit current ring, query, hash, path/subgraph, vector, fingerprint, batch, Python, oracle, and ChEMBL facilities and append a concrete reuse and duplication gap report to `dev/gap_reports/rdkit_layered_fingerprint_source_inventory.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Port complete `_complexQueryHelper`, atom/bond `isComplexQuery`, `isAtomAromatic`, `queryIsBondInRing`, and `queryBondMinRingSize` behavior through narrow typed helpers with line-level markers and exhaustive local tests.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run the exact focused Rust query-complexity, query-aromaticity, and ring-accessor tests added in Step 7 with `op-contracts-strict`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Consolidate the pinned 32-bit `gboost::hash_range` beside the existing sole `hash_combine` with exact empty, order, overflow, and non-native-width tests.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Run the exact focused Rust Boost hash tests added in Step 11 with `op-contracts-strict`.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Consolidate the existing RDKFingerprint linear-path and branched-subgraph enumeration into one reusable crate-private interface with exact absent/empty/duplicated/reordered roots and aggregation-order tests.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Run the exact focused Rust shared path-enumeration test target added in Step 15 with `op-contracts-strict`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port Layered preconditions, SSSR preparation, bond cache, query mask, atomic-number cache, and aromaticity cache with exact preparation-state and intermediate tests using the source-backed query helpers.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run the exact focused Rust Layered preparation and query-mask target added in Step 19 with `op-contracts-strict`.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Port complete Layered `0x01` topology and `0x02` bond-order encoders with exact degree, neighbor, canonicalization, suppression, packing, modulo, and overflow tests.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Run the exact focused Rust Layered topology and bond-order encoder target added in Step 23 with `op-contracts-strict`.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Port complete Layered `0x04` atom-type and `0x20` aromaticity encoders with exact endpoint order, suppression, packing, modulo, and query-aware aromaticity tests.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run the exact focused Rust Layered atom-type and aromaticity encoder target added in Step 27 with `op-contracts-strict`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Port complete Layered `0x08` ring-presence and `0x10` minimum-ring-size encoders with exact omission, zero, suppression, SSSR, modulo, and fused-ring tests.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run the exact focused Rust Layered ring encoder target added in Step 31 with `op-contracts-strict`.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port the complete Layered projection block with exact sort, suffixes, 32-bit hash, modulo, mask, collision, and once-per-accepted-path count tests.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run the exact focused Rust Layered projection target added in Step 35 with `op-contracts-strict`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Port complete `LayeredFingerprintMol` orchestration into one project-native Rust parameter/result API with defaults, typed flags, metadata, rooted enumeration, and focused end-to-end tests.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Run the exact focused Rust Layered end-to-end and public API targets added in Step 39 with `op-contracts-strict`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add project-native Python Layered scalar and ordered batch methods delegating to the sole Rust core with type, default, error, count, mask, root, immutability, order, and concurrency tests.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run the exact focused Python Layered API target added in Step 43.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Add pinned Layered profiles, oracle generator, expected-data schema, version assertion, identity manifest, and focused/small/5,000-row registration under the existing testdata workflow.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run Layered golden preparation for every focused, small, and 5,000-row profile and verify all manifests and checksums.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Add exact Layered focused, complete-small, and complete-5,000-row Rust parity tests for every layer, size, path, root, count, mask, and query profile without row filtering.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Run the exact Layered focused, complete-small, and complete-5,000-row parity targets added in Step 51 in release mode with `op-contracts-strict`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Add source-derived Layered query fixtures and exact intermediate/final parity tests for every bond/endpoint complexity mask and query-aware aromaticity branch.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run the exact focused Rust Layered query parity target added in Step 55 with `op-contracts-strict`.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Extend the repository-owned ChEMBL 37 workflow with complete deterministic Layered profiles, exact schemas, sharding, identity, resumability, and zero-mismatch acceptance tests.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run the exact focused workflow tests for the Layered ChEMBL additions from Step 59.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Run every active Layered profile over every mutually parseable ChEMBL 37 row with deterministic parallel execution and preserve the complete external run identity, counts, checksums, and mismatch artifacts.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Add representative Layered benchmarks against pinned RDKit for path-rich, ring-rich, query-bearing, disconnected, and large molecules with allocation and enumeration measurements.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Run the Layered benchmark suite from Step 65 and preserve its environment, raw results, and complexity review outside tracked machine artifacts.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Add Layered composition tests for repeated, interleaved, and concurrent calls with every implemented fingerprint family on shared molecules without mutation, cache leakage, ordering drift, or option cross-talk.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Run the exact Rust and Python Layered composition targets added in Step 69 with `op-contracts-strict` on Rust.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Add static architecture tests rejecting duplicate Layered enumeration, query helpers, encoders, hash implementations, vectors, and scalar/batch cores.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run the exact static Layered architecture target added in Step 73 with `op-contracts-strict`.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Run `cargo fmt --all`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Regenerate `python/cosmolkit.pyi` with the project stub generator.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Add focused stub-surface tests for Layered scalar and batch signatures.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Run the exact Layered stub-surface tests added in Step 85.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Run `.venv/bin/basedpyright python/tests python/examples`.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Run the focused Layered Python test suite.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Audit Layered source-line, branch, option, output, error, allocation, complexity, API, and corpus closure and write `dev/gap_reports/rdkit_layered_fingerprint_full_port_validation.md` with exact evidence and no incomplete success claim.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Update Layered Rust/Python docs and examples with exact layers, paths, roots, masks, counts, upstream experimental metadata, errors, and composition behavior.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Run the repository documentation build, link, and doctest validation entrypoint for Layered pages and examples.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Update Layered feature metadata, inventory, parity scope, README status, and current Unreleased changelog from the completed validation report without modifying historical plans.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Run repository documentation and status consistency checks for the Layered claims.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Re-run every focused and complete Layered parity target after documentation and generated-artifact changes.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Run the final strict workspace release suite after all Layered artifacts are complete.
