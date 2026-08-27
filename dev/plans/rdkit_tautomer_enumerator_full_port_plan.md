# RDKit Tautomer Enumerator Full Source-Port Plan

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

## Objective And Evidence Boundary

- Port the complete pinned RDKit `MolStandardize::TautomerEnumerator` behavior into one Rust implementation, including transform catalogs, enumeration, modified-atom and modified-bond reporting, termination status, callback cancellation, stereo/isotopic-hydrogen handling, default and custom scoring, canonical selection, and the V1 catalog.
- Chemistry oracle and source tree: RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`, vendored under `third_party/rdkit/`.
- The public API will use COSMolKit naming and value semantics; RDKit CamelCase names remain source anchors and parity-oracle labels rather than a duplicate public API family.
- Exact parity covers ordered canonical isomeric SMILES, full emitted molecular state, enumeration status, modified atom and bond sets, score components, canonical selection, option branches, callbacks, error branches, and input immutability.
- Source defaults remain exact: `max_tautomers = 1000`, `max_transforms = 1000`, `remove_sp3_stereo = true`, `remove_bond_stereo = true`, `remove_isotopic_hydrogens = true`, and `reassign_stereo = true`.
- The existing multi-molecule operation machinery is the only topology-mutation route for enumeration; every candidate branch must complete its own operation contract and molecule-invariant validation before it can enter the result.
- Existing canonical SMILES, SMARTS parsing and matching, kekulization, ring perception, sanitization, stereochemistry, CIP, typed atom/bond identifiers, and molecule storage remain the single shared implementations; this plan must fix a shared prerequisite when parity proves it wrong and must not create tautomer-local substitutes.

## Public API Boundary

```text
TautomerTransform
TautomerCatalog
TautomerOptions
TautomerEnumerator
TautomerEnumeration
TautomerEnumerationStatus
TautomerScoreTerm
TautomerScore

Molecule::enumerate_tautomers()
Molecule::enumerate_tautomers_with_options(...)
TautomerEnumerator::enumerate(...)
TautomerEnumerator::pick_canonical(...)
TautomerEnumerator::canonicalize(...)
```

- `TautomerEnumeration` owns molecules in canonical-SMILES order and exposes status, modified atoms, modified bonds, canonical SMILES, indexed/iterated access, and canonical selection without a second enumeration.
- Rust modified-state sets use `AtomId` and `BondId`; Python exposes integer indices while retaining the same ordering and bounds behavior.
- Enumeration returns a rich typed result through `output: multiple`, `result_type`, and `assemble_fn`; the assembler receives only contract-finalized molecules plus non-molecule metadata produced by the operation body.
- `canonicalize()` delegates to the same enumeration and scoring implementation and returns a new molecule.
- RDKit's `canonicalizeInPlace()` assignment behavior is reproduced by a private compatibility test but is not exposed as a second public mutation API; callers use value assignment, preserving the project's operation and naming policies without losing chemistry behavior.
- Rust callbacks and custom scorers are borrowed callable values for one invocation; Python callbacks receive stable public result views and cannot access internal branch handles.

## Normative Call Graph

```text
Rust/Python project API
  -> registered enumerate_tautomers operation (output: multiple)
     -> TautomerEnumerator configuration and TautomerCatalog
     -> canonical isomeric SMILES + property-cache preparation
     -> symmetrize SSSR + canonical kekulization
     -> ordered canonical-SMILES candidate map
     -> each unfinished candidate
        -> each ordered TautomerTransform
           -> SMARTS SubstructMatch
           -> source-order hydrogen/bond/charge edits
           -> exact partial sanitize flags
           -> stereo and isotopic-H lifecycle
           -> canonical-SMILES deduplication
           -> per-branch operation-contract finalization
     -> limit/callback status and source-defined pruning/rekeying
     -> ordered emitted molecules + metadata assembler
  -> canonical selection
     -> ring score + substructure score + hetero-H score
     -> highest score, then lexicographically smallest canonical SMILES
     -> forced clean stereochemistry assignment on the selected copy
```

## Source Function Ledger

| Source unit | Functions or methods that must be represented exactly | Rust ownership |
|---|---|---|
| `TautomerCatalogParams.h` | `TautomerTransform` constructor, copy constructor, assignment operator | `TautomerTransform` value type; one compiled SMARTS query and source-ordered bond/charge edits |
| `TautomerCatalogUtils.cpp` | both anonymous `getTautomer()` overloads, `stringToBondType()`, `stringToCharge()`, all three `readTautomers()` overloads | one catalog parser used by built-ins, text, stream, and file inputs |
| `TautomerCatalogParams.h/.cpp` | default/file/data/copy constructors, `getTransforms()`, `getTransform()`, `toStream()`, `Serialize()` | `TautomerCatalog`; Rust `Clone`, indexed access, deterministic serialization |
| `TautomerCatalogParams.cpp` | `initFromStream()`, `initFromString()` | explicitly documented upstream-under-construction boundary; no fabricated deserializer |
| `TautomerCatalogEntry.h/.cpp` | constructors, `toStream()`, `Serialize()`, `initFromStream()`, `initFromString()` | written reachability note only; the transform is not serialized upstream and this hierarchy entry is absent from the enumerator call path |
| `Tautomer.cpp` | `SubstructTerm::SubstructTerm()`, `getDefaultTautomerScoreSubstructs()` | `TautomerScoreTerm` and the sole static 12-term default table |
| `Tautomer.cpp/.h` | `scoreRings()`, `scoreSubstructs()`, `scoreHeteroHs()`, inline `scoreTautomer()` | component and aggregate `TautomerScore` functions |
| `Tautomer.h` | both `Tautomer` constructors and internal counters/done flag | private candidate record keyed by canonical isomeric SMILES |
| `Tautomer.h` | `TautomerEnumeratorResult` default/copy constructors, iterator, `begin()`, `end()`, `size()`, `empty()`, `at()`, `operator[]`, `modifiedAtoms()`, `modifiedBonds()`, `status()`, `tautomers()`, `operator()`, `smiles()`, `smilesTautomerMap()`, `fillTautomersItVec()` | `TautomerEnumeration` with one ordered representation and no duplicated molecule collection |
| `Tautomer.h/.cpp` | catalog/cleanup/copy constructors, assignment, every option setter/getter, callback setter/getter | `TautomerOptions` plus `TautomerEnumerator`; Rust value/borrow semantics replace pointer ownership without changing behavior |
| `Tautomer.cpp` | `setTautomerStereoAndIsoHs()` | one private stereo/isotopic-H transition function using shared ring/stereo/CIP components |
| `Tautomer.cpp` | deprecated `enumerate(mol, modifiedAtoms, modifiedBonds)` | private compatibility adapter used only to prove correspondence with the rich result |
| `Tautomer.cpp` | `enumerate(const ROMol&)` | the sole enumerator core driven through `MultiMoleculeOpParts` |
| `Tautomer.cpp/.h` | result overload and iterable template overload of `pickCanonical()` | result-aware and iterable public Rust methods sharing one score/tie-break implementation |
| `Tautomer.cpp` | `canonicalize()`, `canonicalizeInPlace()` | public value-returning canonicalization plus private in-place correspondence test |
| `Tautomer.h` | `tautomerEnumeratorFromParams()`, `getV1TautomerEnumerator()` | ordinary constructors `from_options()` and `v1()` with no duplicate core |
| `Wrap/Tautomer.cpp` | result sequence/property wrappers, callback bridge, constructors, enumerate/canonical/scorer bridges, score functions, score term | idiomatic PyO3 classes and snake_case methods over the Rust core |

## Exact Enumeration State Machine

1. Produce the input's canonical isomeric SMILES, conditionally update its property cache with strict checking disabled, symmetrize SSSR, and canonical-kekulize a copy.
2. Insert the original tautomer and its kekulized form into an ordered map keyed by canonical isomeric SMILES.
3. Visit unfinished candidates in map order and transforms in catalog order; match transforms with ordinary `SubstructMatch` semantics.
4. Check tautomer/transform limits and invoke the callback in the exact source order before applying each accepted transform attempt.
5. Clone the kekulized candidate, move one explicit hydrogen from the first match atom to the last, set source-defined `noImplicit`, apply explicit bond lists or alternating single/double bond changes, and apply source charge deltas.
6. Sanitize with exactly `KEKULIZE | SETAROMATICITY | SETCONJUGATION | SETHYBRIDIZATION | ADJUSTHS`; catch only the source-caught kekulization failure and continue.
7. Record transform endpoints and every changed bond, including bonds additionally changed by sanitization, then apply the exact stereo/isotopic-H transition.
8. Deduplicate by canonical isomeric SMILES, canonical-kekulize accepted products, and retain source counters and done flags.
9. When the global modified atom/bond sets grow, reapply stereo/isotopic-H handling, remove old keys, rekey changed candidates, and merge duplicates in source order.
10. Preserve source status correction, including converting `MaxTautomersReached` to `Completed` when duplicate elimination leaves the final size below the cap.

## Exact Scoring Contract

- Each aromatic ring contributes `100`; an all-carbon aromatic ring contributes an additional `150`.
- The sole built-in table contains the pinned 12 SMARTS/name/weight triples in source order.
- Hydrogens attached to P, S, Se, or Te contribute `-1` each through the source total-hydrogen semantics.
- Canonical selection chooses the greatest signed score; ties choose the lexicographically smallest canonical SMILES; the returned copy receives clean forced stereochemistry assignment.
- Custom scorers replace only the score function, never the canonical-SMILES tie-break or final stereochemistry assignment.

## Existing Infrastructure Integration And Guardrails

- `TautomerTransform` stores the existing compiled SMARTS/query representation and calls the existing matcher; it does not include a local SMARTS parser, VF2 implementation, or reaction/SMIRKS engine.
- Candidate keys use the existing canonical isomeric SMILES writer; no tautomer-local canonicalization or hash key is permitted.
- Ring, kekulization, sanitization, total-hydrogen, stereo, CIP, atom-property, and bond-property transitions delegate to their existing source-backed implementations.
- Catalog constants are generated once from the pinned 37 current and 36 V1 source transforms and are consumed by the same parser path as custom data.
- Repeated matching must reuse compiled transform queries; no per-candidate SMARTS recompilation is allowed.
- No heuristic transform filtering, tautomer pruning, chemistry-specific cap, special-case molecule list, or corpus-conditioned fallback is permitted.
- If a shared prerequisite diverges, the fix belongs in that shared implementation with a local regression and cross-feature non-regression evidence.
- The operation registry may add a narrowly guarded source-defined CIP transition only if strict validation proves `Preserve`, `ClearComputed`, and `Assign` cannot express the option-dependent source behavior; this permission must be scoped to the tautomer operation and must not grant raw storage access.

## Corpus And Branch Matrix

- Focused source cases: every relevant branch in `testTautomer.cpp`, `catch_tests.cpp`, and `Wrap/Tautomer.cpp`, including invalid transforms, charged systems, isotopes, E/Z, tetrahedral centers, ring bonds, limits, callback cancellation, custom scoring, and atom-order permutations.
- Catalog coverage: every one of the 37 current and 36 V1 transforms, including parsed bond and charge vectors and at least one positive match where the source suite provides one.
- Canonical corpus: `third_party/rdkit/rdkit/Chem/MolStandardize/test_data/1kPCS_tautomer.csv.gz` and `100kPCS_tautomer.csv.gz` with exact source endpoint expectations.
- Project corpora: focused project fixtures, the committed 5,000-molecule corpus, and the reproducible ChEMBL 37 pipeline.
- Every corpus row compares parse outcome, ordered tautomer count and strings, status, modified sets, full molecule state, canonical result, score components, all option profiles, and deterministic repeated/threaded execution where applicable.
- The ChEMBL run records exact row/profile/branch/observation counts and zero mismatches before support status is promoted; a small retained replay cannot substitute for the complete run.

## Deliberate Non-Goals

- `Code/GraphMol/TautomerQuery` is a separate query-construction and serialization feature and is not a prerequisite for enumerator correctness; create a follow-up plan only after this plan completes.
- Full `MolStandardize` cleanup, normalization, reionization, fragment-parent, and metal-disconnection pipelines are outside this plan.
- A generic reaction/SMIRKS execution engine is outside this call graph and must not be introduced solely for these ordered atom/bond/charge transforms.
- Broken or upstream-under-construction catalog deserialization is not converted into plausible functionality.
- Deprecated RDKit names and pointer-ownership APIs are not exposed as a second production surface.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit every pinned tautomer source function, source test branch, current COSMolKit prerequisite, and target API mapping and write `dev/gap_reports/rdkit_tautomer_enumerator_source_inventory.md`.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Implement `TautomerTransform` construction, cloning, assignment-equivalent value semantics, compiled-query ownership, and source-ordered bond/charge edit storage as the sole transform representation.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add focused tests for `TautomerTransform` construction, deep clone independence, compiled-query reuse, empty edit vectors, and mismatched transform-shape errors.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::transform`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port `stringToBondType()` and `stringToCharge()` with exact token, whitespace, empty-field, sign, invalid-token, and source integer-conversion behavior.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Add focused catalog-token tests for every source bond symbol, charge form, empty form, whitespace branch, invalid token, and numeric boundary.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::catalog_tokens`.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Port both anonymous `getTautomer()` overloads as one line-to-transform construction path with exact tab-column, comment, default-field, SMARTS parse, bond-vector, and charge-vector behavior.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Add focused transform-line tests for every column count, comment/blank case, omitted optional column, malformed SMARTS, explicit bonds, explicit charges, and source error branch.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::transform_lines`.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Port all three `readTautomers()` overloads into one shared reader with exact file-open, stream-count, comment skipping, order preservation, short-input, and in-memory-data semantics.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Add focused reader tests for file, stream, bounded stream, embedded definitions, zero limit, early EOF, comments, blank lines, deterministic order, and structured I/O failures.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::catalog_reader`.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Port the current 37-transform and V1 36-transform built-in tables through generated source-owned constants that feed the sole catalog construction path.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Add exact built-in catalog tests for every transform name, SMARTS, bond vector, charge vector, compiled query, source order, current count, and V1 count.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::builtin_catalogs`.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port `TautomerCatalogParams` default/file/data/copy constructors, `getTransforms()`, `getTransform()`, `toStream()`, and `Serialize()` into `TautomerCatalog` without a second transform store.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Add focused catalog-object tests for all constructors, indexing bounds, clone independence, exact count-only serialization, and the structured unsupported deserialization boundary.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::catalog_object`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Document `TautomerCatalogParams::{initFromStream,initFromString}` and `TautomerCatalogEntry` serialization as unreachable or upstream-under-construction source boundaries in the tautomer source inventory.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Port `SubstructTerm::SubstructTerm()` and `getDefaultTautomerScoreSubstructs()` into `TautomerScoreTerm` and one lazily compiled exact 12-term table.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Add exact scoring-term tests for all names, SMARTS, signed weights, equality semantics, invalid SMARTS handling, source order, and one-time compilation.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::score_terms`.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Port `scoreRings()` with exact SymmSSSR fallback, aromatic-bond and all-carbon bitsets, bond-ring traversal, and `100 + 150` accumulation.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Add focused ring-score tests for no rings, aromatic heterorings, all-carbon aromatics, fused rings, nonaromatic rings, precomputed/non-precomputed SymmSSSR, and source immutability.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::ring_score`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Port `scoreSubstructs()` with exact all-match counting, invalid-matcher skipping, signed multiplication, custom-term order, and shared matcher semantics.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Add focused substructure-score tests for zero/one/multiple/overlapping matches, every built-in term, custom positive/negative terms, invalid matchers, and integer accumulation boundaries.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::substructure_score`.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Port `scoreHeteroHs()` and inline `scoreTautomer()` with exact P/S/Se/Te total-hydrogen counting and signed component aggregation.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Add focused hetero-H and aggregate-score tests for implicit/explicit/isotopic hydrogens, charged atoms, each penalized element, excluded elements, and exact component sums.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::aggregate_score`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Port both private `Tautomer` constructors and all candidate counters/done-state transitions into one canonical-SMILES-keyed candidate record.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Add focused candidate-record tests for default state, explicit state, clone independence, counter transitions, and deterministic ordered-map replacement.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::candidate_record`.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Port every `TautomerEnumeratorResult` constructor, iterator, accessor, indexed lookup, molecule/smiles projection, map view, modified-set view, status view, and iterator-vector behavior into one `TautomerEnumeration` representation.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Add focused result tests for empty/single/multiple results, canonical-SMILES ordering, forward/backward iteration, clone independence, index bounds, projections, typed modified sets, and every status.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::enumeration_result`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Port the catalog/cleanup/copy constructors, assignment semantics, six option setter/getter pairs, and callback setter/getter into `TautomerOptions` and `TautomerEnumerator` with exact defaults and ownership behavior.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Add focused enumerator-configuration tests for defaults, every option boundary, zero limits, clone independence, catalog sharing/value behavior, callback replacement, and callback lifetime.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::enumerator_configuration`.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Audit tautomer option-dependent CIP transitions against strict operation validation and write the exact required policy mapping into `dev/gap_reports/rdkit_tautomer_enumerator_source_inventory.md`.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Modify the operation contract only if required by Step 85 to add a tautomer-scoped source-defined CIP transition with compile-time registry guardrails and no new raw-storage capability.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Add strict operation tests for every allowed tautomer CIP transition and for rejection by every non-tautomer operation.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict operations::ops::tests::tautomer_cip`.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Port `setTautomerStereoAndIsoHs()` with every modified-atom chiral/CIP branch, isotopic-H branch, modified-bond stereo/direction branch, ring fallback, forced reassignment branch, and `_StereochemDone` branch.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Add focused stereo/isotopic-H tests for every option combination, SP2/SP3 atoms, present/absent CIP properties, explicit isotopic H, ring/non-ring double bonds, E/Z/ANY/NONE, neighboring directions, fast-ring fallback, and source immutability.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::stereo_and_isotopic_hydrogens`.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Port the `enumerate()` input-preparation and initial-candidate source block as one complete helper covering canonical isomeric SMILES, non-strict property-cache update, SymmSSSR, canonical kekulization, ordered insertion, and initialization failures.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Add focused initialization tests for sanitized/unsanitized inputs, missing/present caches, aromatic systems, empty/disconnected molecules, canonical atom-order permutations, initial kekulization errors, and input immutability.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::enumeration_initialization`.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Port the `enumerate()` single-match transformation source block as one complete helper covering hydrogen transfer, `noImplicit`, explicit/alternating bond edits, charge deltas, exact sanitize flags, caught/propagated errors, modified sets, stereo handling, dedup key, and canonical kekulization.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Add focused transformation tests for every built-in edit shape, endpoint order, explicit/implicit/isotopic hydrogen, alternating and explicit bonds, charge deltas, sanitization-changed bonds, duplicate products, caught kekulization failure, propagated non-kekulization failure, and per-branch invariants.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::single_transform_application`.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Port the `enumerate()` ordered expansion loop as one complete helper covering unfinished-candidate traversal, transform/match order, source counter increments, both limits, callback timing, insertion, done flags, and termination statuses.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Add focused expansion-loop tests for map insertion during traversal, multiple matches, duplicate matches, zero/exact/over limits, transform count semantics, callback observation/cancellation order, done flags, and deterministic termination.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::enumeration_expansion`.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Port the `enumerate()` post-expansion pruning/rekeying source block as one complete helper covering modified-set growth, stereo reapplication, key replacement, duplicate merging, count updates, final status correction, and iterator-order materialization.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Add focused pruning tests for unchanged/changed stereo keys, duplicate collapse, ordered rekeying, modified-set growth across rounds, `MaxTautomersReached` correction, count consistency, and final iterator order.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::enumeration_pruning`.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Port `TautomerEnumerator::enumerate(const ROMol&)` by composing only the completed source-block helpers through `MultiMoleculeOpParts` and assembling one `TautomerEnumeration`.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Add operation-level enumeration tests for registry metadata, emitted-branch validation, source/child derivation, foreign-handle rejection, output order, rich-result assembly, input COW/value semantics, coordinates/properties, and strict invariants.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test tautomer_enumeration_operation`.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Port the deprecated vector-returning `enumerate()` overload as a private compatibility adapter over `TautomerEnumeration` with exact optional modified-bitset output behavior.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Add compatibility tests proving the private deprecated adapter equals every rich-result molecule and modified set without exposing a duplicate production algorithm.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::tautomer::tests::deprecated_enumerate_correspondence`.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Port both `pickCanonical()` overloads with exact single-item optimization, signed score maximum, canonical-SMILES tie-break, result-map optimization, iterable behavior, custom scorer, returned-copy semantics, and final clean forced stereochemistry assignment.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Add focused canonical-picking tests for empty/single/multiple inputs, positive/negative/equal scores, lexical ties, result/iterable equivalence, custom scorers, source immutability, and final stereochemistry.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::canonical_selection`.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Port `canonicalize()`, private `canonicalizeInPlace()` correspondence behavior, `tautomerEnumeratorFromParams()`, and `getV1TautomerEnumerator()` as delegates of the sole enumerator/catalog/scoring implementation.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Add focused canonicalization and factory tests for default/current/V1 catalogs, custom options, custom scorer, atom-order-independent endpoint selection, private in-place correspondence, and absence of a duplicate public mutation route.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict chemistry::tautomer::tests::canonicalization_and_factories`.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Implement the project-native Rust exports, molecule methods, error types, support descriptor, operation registry entries, invariant matrix entries, and parity matrix entries as thin surfaces over the sole tautomer core.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Add public Rust API tests for discoverability, defaults, options, result ergonomics, custom catalogs, callbacks, custom scoring, structured errors, composition, thread safety, and absence of RDKit-style duplicate names.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test tautomer_public_api`.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Implement pinned-RDKit reference generators and committed manifest workflows for focused, catalog, PCS, 5,000-row, and ChEMBL 37 tautomer evidence under `tools/testdata/rdkit/` and `dev/tools/chembl_parity/`.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Generate the committed focused and catalog oracle fixtures with source revision, RDKit version, generator command, schema, checksums, profile definitions, and row provenance.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Add exact focused parity tests covering every reachable `testTautomer.cpp`, `catch_tests.cpp`, and wrapper branch plus every current/V1 transform.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_tautomer_focused_parity`.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Add exact `1kPCS_tautomer.csv.gz` and `100kPCS_tautomer.csv.gz` canonical endpoint, input-tautomer invariance, atom-order permutation, score, and full-enumeration tests.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_tautomer_pcs_parity`.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Add exact project-corpus tests for focused, 5,000-row, and composition profiles comparing parse outcome, ordered outputs, full molecular state, status, modified sets, score components, canonical result, options, repeated calls, and parallel batches.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_tautomer_corpus_parity`.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [ ]: Run the complete ChEMBL 37 tautomer matrix with source defaults and every declared option/composition profile and write its exact aggregate evidence artifact.
Step 170 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [ ]: Fix every ChEMBL 37 mismatch at the corresponding pinned source function or shared prerequisite and add each new minimal failing molecule to the owning focused regression fixture.
Step 172 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [ ]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_tautomer_focused_parity`.
Step 174 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [ ]: Rerun the complete ChEMBL 37 tautomer matrix after all source-level fixes and replace the evidence artifact only when every declared comparison is exact.
Step 176 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [ ]: Audit repeated SMARTS matching, candidate cloning, ordered-map operations, canonical SMILES, sanitization, stereo reassignment, and branch finalization against RDKit allocation and complexity shape and write the findings into the source inventory.
Step 178 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [ ]: Fix every proven tautomer hot-path complexity or allocation regression by reproducing the source data flow without heuristic pruning or behavior changes.
Step 180 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 181 [ ]: Add deterministic performance fixtures that compare source and Rust scaling across transform count, tautomer count, atom count, ring complexity, and repeated SMARTS matches without asserting unstable wall-clock thresholds.
Step 182 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 183 [ ]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test tautomer_complexity_regression`.
Step 184 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 185 [ ]: Implement the PyO3 enumerator, options, result sequence, status, score-term, callback, custom-scorer, molecule-method, and exception surfaces as thin snake_case wrappers over the Rust API.
Step 186 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 187 [ ]: Regenerate `python/cosmolkit.pyi` using the project stub generator.
Step 188 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 189 [ ]: Add Python tests for default/custom enumeration, result iteration/indexing, modified sets, statuses, catalog data, callback cancellation, custom scoring, canonical selection, exception mapping, input immutability, and cross-feature composition.
Step 190 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 191 [ ]: Run `.venv/bin/pytest python/tests -k tautomer`.
Step 192 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 193 [ ]: Update current public Rust/Python documentation, examples, support matrices, `README.md`, `VALIDATION.md`, and the current release section with the completed tautomer API and exact final evidence without editing historical plans.
Step 194 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 195 [ ]: Add documentation and example smoke tests for parse-to-enumerate-to-canonical-to-SMILES/SDF/InChI/fingerprint/descriptor composition in Rust and Python.
Step 196 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 197 [ ]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test tautomer_documented_workflows`.
Step 198 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 199 [ ]: Run `.venv/bin/pytest python/tests -k tautomer_documented_workflows`.
Step 200 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 201 [ ]: Run `cargo fmt --all -- --check`.
Step 202 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 203 [ ]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 204 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 205 [ ]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 206 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 207 [ ]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.
Step 208 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 209 [ ]: Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml`.
Step 210 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 211 [ ]: Run `.venv/bin/pytest`.
Step 212 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 213 [ ]: Run `.venv/bin/basedpyright python/tests python/examples`.
Step 214 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 215 [ ]: Run `.venv/bin/python -m sphinx -W --keep-going -b html python/docs/source python/docs/build/html`.
Step 216 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 217 [ ]: Audit every tautomer source marker, function-ledger row, operation-matrix row, corpus manifest, public symbol, Python stub, documentation statement, example, and final evidence count and write the completion result into `dev/gap_reports/rdkit_tautomer_enumerator_source_inventory.md`.
