# RDKit Topological Torsion Fingerprint Full Source-Port Plan

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
- Before every task that edits `crates/cosmolkit-core/`, read `dev/README.md` in the step immediately before the required policy/source reading step.
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

## Source Pin

- Chemistry oracle and source: RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- The source pin is the repository-wide reference recorded in `testdata/reference/rdkit.json`; moving to another RDKit revision requires regenerating this plan's source inventory and parity data.
- Source URLs in the inventory must use the pinned revision rather than a moving branch.
- The source license is RDKit BSD-3-Clause; copied source lines remain review comments governed by `dev/source_reproduction_protocol.md`.

## Scope

The port covers the complete pinned-RDKit Topological Torsion fingerprint family reached through both modern and legacy public entry points:

- modern `TopologicalTorsion::getTopologicalTorsionGenerator<std::uint64_t>()` construction, mutable options, four scalar vector forms, four multi-molecule vector forms, custom atom invariants, `fromAtoms`, `ignoreAtoms`, chirality, `onlyShortestPaths`, `numBitsPerFeature`, count simulation, all supported `AdditionalOutput`, information strings, and JSON roundtrips;
- legacy `AtomPairs::getTopologicalTorsionFingerprint()`, `getHashedTopologicalTorsionFingerprint()`, and `getHashedTopologicalTorsionFingerprintAsBitVect()` behavior, including the un-hashed vector-size compatibility bug and the distinct non-four-bit count-threshold rule;
- Python generator construction and scalar/bulk calls, legacy descriptor wrappers, and the public `rdkit.Chem.AtomPairs.Torsions` score/explanation/id helpers plus `AtomPairs.Utils.ExplainAtomCode`, mapped to COSMolKit-native snake-case names without duplicating the chemistry kernel;
- Rust scalar and batch convenience APIs, feature metadata, docs, examples, generated stubs, deterministic oracle generation, focused branch parity, and maintained-corpus exact parity.

The port does not include the Atom Pair fingerprint family merely because Topological Torsion reuses `AtomPairAtomInvGenerator`. Only the shared atom-code and atom-invariant source unit reached by Topological Torsion is in scope. The operation registry is not involved because fingerprint generation is read-only and does not mutate public `Molecule` values.

## Single-Core Architecture Contract

- `TopologicalFingerprintParams` and `Molecule::topological_fingerprint()` remain the existing RDKit path/subgraph fingerprint (`RDKFingerprintMol`) API; Topological Torsion receives distinct `TopologicalTorsion*` names.
- `AdditionalOutput`, `FingerprintArguments`, `FingerprintFuncArguments`, `SparseCountFingerprint`, `SparseBitFingerprint`, `Fingerprint`, the Boost-compatible `hash_combine`, and the RDKit fingerprint Mersenne Twister remain single shared implementations.
- The source-backed `FingerprintGenerator.cpp` output driver currently embedded in `MorganFingerprintGenerator` must be extracted into one private shared generator core before Topological Torsion is connected; that extraction must also reproduce the exact missing-`_StereochemDone` private-copy preparation and source-object selection flow without mutating the caller's molecule or adding a fingerprint-local stereochemistry implementation.
- Morgan and Topological Torsion must delegate sparse-count, sparse-bit, folded-count, explicit-bit, count-simulation, random extra-bit, and `AdditionalOutput` propagation to that same shared core.
- The unused historical `build_fingerprint()`/`fold_invariant()`/`atom_is_excluded()` Morgan branch must be removed before the torsion port so it cannot become an accidental second output path.
- The existing Boost-compatible `hash_combine()` is reused directly by `getTopologicalTorsionHash()`; no torsion-local hash implementation is allowed.
- The existing source-backed `MolOps::getDistanceMat`/`FloydWarshall` logic currently private to distance geometry must be moved to one shared internal graph-distance helper and reused by `onlyShortestPaths`; a second Floyd-Warshall implementation is forbidden.
- Existing RDKFingerprint path/subgraph enumeration must not be reused as a shortcut because it ports different upstream functions and returns bond paths with different branching and rooting semantics.
- Legacy and modern public adapters may remain distinct where their upstream size, validation, warning, or threshold behavior differs, but they must terminate in the same atom invariant, path enumeration, torsion code/hash, environment, and vector-assembly kernels.
- No production path may call RDKit, a Python oracle, fixture-specific lookup table, or historical approximate fingerprint code.

## Complete Calling Path

```text
Modern Python/Rust generator entry
└─ GetTopologicalTorsionGenerator / get_topological_torsion_generator
   └─ TopologicalTorsionArguments + AtomPairAtomInvGenerator(correction=true)
      └─ shared FingerprintGenerator core
         ├─ getSparseCountFingerprint -> getFingerprintHelper(fpSize=0)
         ├─ getSparseFingerprint -> getFingerprintHelper(effectiveSize) -> count simulation
         ├─ getCountFingerprint -> getFingerprintHelper(fpSize)
         └─ getFingerprint -> getFingerprintHelper(effectiveSize) -> count simulation
            ├─ optional assignStereochemistry/CIP prerequisite
            ├─ AtomPairAtomInvGenerator::getAtomInvariants
            │  └─ getAtomCode -> numPiElectrons + CIP label lookup
            └─ TopologicalTorsionEnvGenerator::getEnvironments
               ├─ findAllPathsOfLengthN
               │  └─ findAllPathsOfLengthsMtoN
               │     ├─ shared MolOps distance matrix when onlyShortestPaths=true
               │     └─ pathFinderHelper -> extendPaths -> bond-set deduplication
               ├─ endpoint fromAtoms filter + whole-path ignoreAtoms filter
               ├─ repeated-atom/ring-closure acceptance rule
               ├─ endpoint/internal branch correction
               ├─ getTopologicalTorsionCode when fpSize=0
               └─ getTopologicalTorsionHash -> shared hash_combine when fpSize!=0

Legacy un-hashed entry
└─ getTopologicalTorsionFingerprint
   └─ same invariant/path/code kernels -> legacy result size (2^width - 1)

Legacy hashed-count entry
└─ getHashedTopologicalTorsionFingerprint
   └─ TorsionFpCalc -> modern generator getCountFingerprint

Legacy hashed-bit entry
└─ getHashedTopologicalTorsionFingerprintAsBitVect
   └─ TorsionFpCalc(blockLength) -> legacy nBitsPerEntry threshold expansion
```

## Existing COSMolKit Integration Ledger

| Responsibility | Current implementation | Mainline decision |
|---|---|---|
| Optional provenance containers | `fingerprint.rs::AdditionalOutput` | Reuse; correct repeated-use behavior against pinned `reinitAdditionalOutput`, including the source's `atomsPerBit` handling. |
| Common generator arguments and per-call arguments | `FingerprintArguments`, `FingerprintFuncArguments` | Reuse; extend only with torsion-specific typed arguments outside the common structure. |
| Sparse count, sparse bit, explicit bit vectors | `SparseCountFingerprint`, `SparseBitFingerprint`, `Fingerprint` | Reuse unchanged unless source-width tests expose a concrete gap. |
| Generic count simulation and `AdditionalOutput` duplication | Methods on `MorganFingerprintGenerator`, `duplicate_additional_output_bit`, `setup_temp_additional_output` | Extract once into a private shared generator core; retain one implementation. |
| Shared chirality preparation | The active generator helper currently reports an unsupported option when chirality is requested without pre-existing `_StereochemDone` state | Complete the pinned `getFingerprintHelper` private-copy `MolOps::assignStereochemistry` branch with exact `mol`/`lmol` selection and delegate to the repository's shared stereo machinery; do not add a torsion-local stereo implementation or silently change which object supplies invariants. |
| Extra-bit RNG | `RdkitFingerprintMtRng` | Reuse through the shared generator core; test torsion `numBitsPerFeature > 1`. |
| Boost 32-bit hash combine | `fingerprint.rs::hash_combine` | Reuse directly; do not add a second formula. |
| Topological distance matrix | Private source-backed implementation in `chemistry/distgeom.rs` | Move to a shared private graph helper and reuse from both callers. |
| RDKFingerprint path enumeration | `enumerate_rdkit_fp_paths` and helpers | Keep separate because it corresponds to `RDKitFPUtils::enumerateAllPaths`, not Topological Torsion's `findAllPathsOfLengthN`. |
| Historical Morgan output construction | Unused `build_fingerprint`, `fold_invariant`, `atom_is_excluded` | Delete before introducing torsion so only the active generator core remains. |
| Existing `topological_fingerprint` public family | RDKit `RDKFingerprintMol` | Preserve without renaming or semantic reuse; add a distinct Topological Torsion family. |

## Per-Function Source Ledger

| Source file and pinned lines | Function or source unit | Planned artifact |
|---|---|---|
| `Code/GraphMol/Fingerprints/FingerprintUtil.h:28-41` | AtomPairs bit widths, atom-number buckets, and path constants | One shared constant set used by atom explanation, atom coding, legacy sizing, and torsion coding. |
| `Code/GraphMol/Atom.cpp:934-953` | `numPiElectrons` | Shared source-backed atom helper reused by torsion atom coding. |
| `Code/GraphMol/Chirality.cpp:2889-2905` and `Code/GraphMol/MolOps.h:1185-1187` | fingerprint-reached default `MolOps::assignStereochemistry` call | Shared private-copy prerequisite with exact early return/default flags and an explicitly audited observable boundary; broader unrelated stereo APIs remain outside this plan. |
| `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:45-97` | `AtomPairs::getAtomCode` | Shared AtomPairs invariant helper with exact branch/pi/type/chirality packing. |
| `Code/GraphMol/Fingerprints/AtomPairGenerator.h:22-47` | `AtomPairAtomInvGenerator` contract | Shared constructor defaults, polymorphic invariant contract, clone, and metadata surface retained as one typed Rust unit. |
| `Code/GraphMol/Fingerprints/AtomPairGenerator.cpp:26-66` | `AtomPairAtomInvGenerator::{constructor,getAtomInvariants,infoString,toJSON,fromJSON,clone}` | One reusable invariant generator source unit with topological-torsion correction. |
| `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:109-138` | `getTopologicalTorsionCode` | Direction-canonical un-hashed 64-bit torsion encoder. |
| `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:140-167` | `getTopologicalTorsionHash` | Direction-canonical 32-bit hash using the existing shared `hash_combine`. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:189-242` | `Subgraphs::extendPaths` | Ordered atom-path extension with final ring closure behavior and shortest-path pruning. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:244-282` | `Subgraphs::pathFinderHelper` | Ordered path seeding and length-range accumulation. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:442-549` | `findAllPathsOfLengthsMtoN` | Atom adjacency creation, optional distance pruning, atom-to-bond conversion, and bond-set deduplication. |
| `Code/GraphMol/Subgraphs/Subgraphs.cpp:550-555` | `findAllPathsOfLengthN` | Single-length adapter used by torsion generation. |
| `Code/GraphMol/Subgraphs/Subgraphs.h:39-48,107-130` | ordered path types and public path signatures/defaults | Internal Rust ordered-container and argument contract preserving source order and Topological Torsion call defaults. |
| `Code/GraphMol/Matrices.cpp:39-100,166-230` | `FloydWarshall` and `MolOps::getDistanceMat` reached by shortest-path enumeration | Existing source-backed graph-distance kernel moved once and shared without a second matrix implementation. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.h:27-82,292-425` | `AdditionalOutput`, `FingerprintFuncArguments`, generator ownership/options, scalar overloads, and JSON declarations | Existing Rust shared data/ownership surface retained and reconciled with the complete source driver. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:36-87` | `FingerprintArguments::{constructor,commonArgumentsString,toJSON,fromJSON}` | Existing shared implementation retained and regression-checked. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:90-173` | `FingerprintGenerator::{constructor,destructor,infoString}` and `reinitAdditionalOutput` | Shared Rust ownership core and exact output reset semantics. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:182-321` | generator JSON serialization/deserialization | Shared typed JSON dispatch extended with Topological Torsion variants. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:323-435` | `getFingerprintHelper` | One shared environment-to-count-vector driver for Morgan and Topological Torsion. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:438-500` | `duplicateAdditionalOutputBit`, `setupTempAdditionalOutput` | Existing shared helpers retained as the sole count-simulation provenance implementation. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:503-652` | four scalar output methods | One shared sparse/count/bit assembly implementation. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:654-756` | `mtgetFingerprints` and four generator bulk methods | Ordered, deterministic shared batch output implementation. |
| `Code/GraphMol/Fingerprints/FingerprintGenerator.cpp:826-989` | `get*FP` and `get*FPBulk` dispatch | Topological Torsion enum dispatch without family-specific output duplication. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:25-65` | `TopologicalTorsionArguments::{constructor,infoString,toJSON,fromJSON}` | Typed torsion arguments and options. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:33-45` | `TopologicalTorsionEnvGenerator::getResultSize` | Exact output width with checked treatment of source undefined shifts. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:67-99` | `TopologicalTorsionAtomEnv::{updateAdditionalOutput,getBitId}` | Atom-path provenance and stored bit-id behavior. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:101-194` | `TopologicalTorsionEnvGenerator::getEnvironments` | Complete filtering, cycle, correction, code/hash, and source-order environment generation. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:196-210` | environment generator `infoString`, `toJSON`, `fromJSON` | Generator metadata serialization. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.cpp:212-252` | both `getTopologicalTorsionGenerator` overloads | One Rust factory path with safe ownership and 64-bit output restriction. |
| `Code/GraphMol/Fingerprints/TopologicalTorsionGenerator.h:21-136` | argument defaults, atom-environment storage, virtual generator contract, factory defaults, and supported output declaration | Typed Rust surface whose behavior-bearing inline constructor and defaults terminate in the same implementation units listed above. |
| `Code/GraphMol/Fingerprints/AtomPairs.h:197-271` | three deprecated legacy declarations and defaults | Distinct compatibility API defaults routed to the shared chemistry kernels. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:159-265` | `getTopologicalTorsionFingerprint` | Legacy un-hashed adapter with `2^width - 1` compatibility size and shared chemistry kernels. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:267-297` | `TorsionFpCalc` | Legacy-to-modern hashed-count adapter. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:298-310` | `getHashedTopologicalTorsionFingerprint` | Legacy hashed sparse-count result adapter. |
| `Code/GraphMol/Fingerprints/AtomPairs.cpp:312-347` | `getHashedTopologicalTorsionFingerprintAsBitVect` | Legacy block sizing and count-threshold bit expansion. |
| `Code/GraphMol/Fingerprints/Wrap/TopologicalTorsionWrapper.cpp:21-92` | generator Python construction and options | Python typed constructor/options mapping. |
| `Code/GraphMol/Fingerprints/Wrap/FingerprintGeneratorWrapper.cpp:43-320,327-772` | scalar/bulk Python generator methods, options, `AdditionalOutput`, JSON | Reuse one Python-facing vector/provenance surface across fingerprint families. |
| `Code/GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp:964-985` | `AtomPairsParameters` constants and `GetAtomPairAtomCode` | Read-only Python constant/code surface backed by the same shared AtomPairs constant set and atom-code kernel. |
| `Code/GraphMol/Descriptors/Wrap/rdMolDescriptors.cpp:250-309,1035-1070` | three legacy Python wrappers | COSMolKit-native compatibility functions with exact defaults and structured errors. |
| `rdkit/Chem/AtomPairs/Utils.py:18-107` | `ExplainAtomCode` and its `GetAtomCode` dependency | Pure atom-code explanation helper backed by the same shared atom-code constants and kernel. |
| `rdkit/Chem/AtomPairs/Torsions.py:31-162` | aliases, `pyScorePath`, `ExplainPathScore`, `GetTopologicalTorsionFingerprintAsIds` | Pure helper utilities backed by the same atom-code/torsion-code kernel. |

## Edge and Error Closure

The implementation and oracle profiles must explicitly cover empty molecules, zero and one atom, `torsionAtomCount` zero/one/default/maximum, source-defined and source-undefined shift widths, `fpSize` zero, empty and oversized count bounds, `numBitsPerFeature` zero and greater than one, `nBitsPerEntry` zero/one/four/non-four, `nBits` not divisible by `nBitsPerEntry`, empty-versus-null atom selections, duplicate and out-of-range indices, custom-invariant length and unsigned overflow, disconnected graphs, explicit hydrogens, full-ring closures, non-initial repeated atoms, only-shortest-path pruning, query atoms/bonds, dative bonds, aromatic and hypervalent atoms, high-degree modulo behavior, uncommon atomic numbers, stereochemistry-not-assigned input, legacy/non-legacy stereo perception, CIP R/S/unassigned states, hash collisions, count collisions, repeated reuse of one `AdditionalOutput`, JSON missing/unknown/malformed types, bulk ordering, null/invalid batch records, and deterministic multi-thread execution.

Undefined C++ behavior, unchecked invalid memory access, divide-by-zero, and impossible allocation sizes are not reproduced as Rust panics or fabricated chemistry. The execution audit must classify each such input and map it to a structured `FingerprintError` while preserving every source-defined valid branch exactly.

## Completion Criteria

- Every active source line in the function ledger has an inline verbatim two-axis marker at its corresponding Rust implementation or is explicitly classified as a non-behavioral ownership/lifetime difference.
- Morgan, RDKFingerprint, MACCS, and Avalon retain their documented exact outputs after shared-core consolidation.
- Topological Torsion modern sparse-count, sparse-bit, hashed-count, and explicit-bit vectors match pinned RDKit exactly for every focused profile and maintained-corpus row.
- Legacy un-hashed, hashed-count, and hashed-bit outputs match pinned RDKit exactly, including legacy-only size and threshold semantics.
- `atomToBits`, `atomCounts`, `bitPaths`, and `atomsPerBit` match exactly before and after count simulation, collisions, repeated output-object reuse, and bulk calls.
- `fromAtoms`, `ignoreAtoms`, custom atom invariants, chirality, shortest-path behavior, path ordering, ring closure, JSON, information strings, and all public defaults match the selected source boundary.
- Chirality-enabled generation on molecules lacking prior stereochemistry assignment follows the pinned private-copy preparation and exact `mol`/`lmol` source-selection path, matches oracle invariants and outputs, and leaves the source molecule unchanged.
- Rust and Python expose distinct Topological Torsion names and never alias the existing RDKFingerprint `topological_fingerprint` implementation.
- Only one production implementation exists for common vector storage, Boost hashing, random extra-bit generation, generic fingerprint assembly, count simulation, provenance duplication, and topological distance matrices.
- The oracle generator is deterministic, pinned, checksummed, and production code has no oracle dependency.
- Focused tests, 5,000-row exact parity, strict core tests, strict workspace tests, Python tests, generated stubs, docs, examples, and type checks all pass.
- The final validation report records zero unexplained mismatches, source-line closure, algorithmic-complexity review, allocation review, and deterministic benchmark results.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.

Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 3 [x]: Audit the pinned Topological Torsion modern and legacy call graphs against current COSMolKit and write `dev/gap_reports/rdkit_topological_torsion_fingerprint_source_inventory.md` with exact source ranges, public-option mapping, prerequisite closure, error ledger, line-coverage ledger, and same-source reuse decisions.

Step 4 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 6 [x]: Remove the unreachable historical Morgan `build_fingerprint`, `fold_invariant`, and `atom_is_excluded` branch from `crates/cosmolkit-core/src/properties/fingerprint.rs`.

Step 7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 8 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test fingerprint_interaction_parity` after the historical-branch removal.

Step 9 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 11 [x]: Refactor `FingerprintGenerator.cpp:323-652` from `MorganFingerprintGenerator` into one complete private shared Rust generator core that reproduces the private-copy stereochemistry preparation and exact `mol`/`lmol` selection flow without changing Morgan's public API.

Step 12 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 13 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 14 [x]: Add focused strict tests for the shared generator core covering all four scalar outputs, count simulation, collisions, random extra bits, every `AdditionalOutput` allocation, repeated output reuse, missing-`_StereochemDone` chirality preparation, source-molecule immutability, and unchanged Morgan results.

Step 15 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 16 [x]: Run the exact shared-generator-core strict test target added in Step 14.

Step 17 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 19 [x]: Move the existing source-backed `MolOps::getDistanceMat` and `FloydWarshall` Rust implementation from distance geometry into one shared private graph-distance helper.

Step 20 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 21 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 22 [x]: Add focused strict tests proving the shared graph-distance helper preserves distance-geometry matrices, disconnected-pair infinity, source ordering, and shortest-path distances.

Step 23 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 24 [x]: Run the exact shared graph-distance strict test target added in Step 22.

Step 25 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 27 [x]: Implement pinned `Atom.cpp::numPiElectrons` as the single shared atom pi-electron helper with verbatim two-axis source anchors.

Step 28 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 29 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 30 [x]: Implement pinned `FingerprintUtil.cpp::getAtomCode` with exact degree subtraction, modulo widths, atom-type fallback, pi packing, CIP handling, and unsigned arithmetic.

Step 31 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 33 [x]: Implement the complete pinned `AtomPairAtomInvGenerator` source unit with the Topological Torsion correction and shared atom-code delegation.

Step 34 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 35 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 36 [x]: Add focused strict atom-invariant tests for every element bucket, aromaticity, hybridization, explicit hydrogen, dative bond, branch/pi modulo boundary, custom invariant correction, R/S chirality, missing CIP state, and wrapping subtraction.

Step 37 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 38 [x]: Run the exact Topological Torsion atom-invariant strict test target added in Step 36.

Step 39 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 41 [x]: Implement pinned `FingerprintUtil.cpp::getTopologicalTorsionCode` with direction canonicalization, chirality-dependent shifts, and exact 64-bit behavior.

Step 42 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 43 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 44 [x]: Implement pinned `FingerprintUtil.cpp::getTopologicalTorsionHash` by delegating each canonical path code to the existing shared Boost-compatible `hash_combine`.

Step 45 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 47 [x]: Add focused strict torsion-code tests for forward/reverse equality, palindromes, first-difference ordering, chirality widths, zero codes, maximum valid widths, source-undefined shifts, and Boost hash overflow.

Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 49 [x]: Run the exact torsion-code and torsion-hash strict test target added in Step 47.

Step 50 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 51 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 52 [x]: Implement pinned `Subgraphs.cpp::extendPaths` with source-order adjacency scanning, shortest-distance pruning, and final-only ring closure.

Step 53 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 55 [x]: Implement pinned `Subgraphs.cpp::pathFinderHelper` with rooted/unrooted seeding and exact length-map insertion order.

Step 56 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 57 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 58 [x]: Implement pinned `Subgraphs.cpp::findAllPathsOfLengthsMtoN` with hydrogen handling, shared distance-matrix reuse, atom-to-bond conversion, and bond-bitset deduplication.

Step 59 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 61 [x]: Implement pinned `Subgraphs.cpp::findAllPathsOfLengthN` as the single-length adapter over the shared range implementation.

Step 62 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 63 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 64 [x]: Add focused strict ordered-path tests for chains, branches, rings, fused rings, repeated-atom rejection, explicit hydrogens, disconnected graphs, rooted paths, target-length boundaries, bond-set duplicates, and `onlyShortestPaths`.

Step 65 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 66 [x]: Run the exact Topological Torsion path-enumeration strict test target added in Step 64.

Step 67 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 69 [x]: Implement the complete pinned `TopologicalTorsionArguments` source unit with common arguments, torsion atom count, shortest-path option, defaults, information string, and JSON behavior.

Step 70 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 71 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 72 [x]: Implement pinned `TopologicalTorsionEnvGenerator::getResultSize` with exact valid shifts and structured errors for source-undefined widths.

Step 73 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 75 [x]: Add focused strict argument tests for all defaults, option mutations, JSON partial updates, malformed JSON, unknown types, information strings, result sizes, and invalid width mappings.

Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 77 [x]: Run the exact Topological Torsion argument and JSON strict test target added in Step 75.

Step 78 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 79 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 80 [x]: Implement the complete pinned `TopologicalTorsionAtomEnv` source unit with stored bit identifiers and all four supported provenance outputs.

Step 81 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 83 [x]: Implement pinned `TopologicalTorsionEnvGenerator::getEnvironments` with exact source ordering, selections, cycle checks, invariant correction, and code-versus-hash dispatch.

Step 84 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 85 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 86 [x]: Implement the pinned Topological Torsion environment generator information-string and JSON source unit.

Step 87 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 89 [x]: Add focused strict environment tests for null/empty selections, endpoint-only roots, ignored internal atoms, custom invariants, full-ring closures, illegal repeated atoms, hashed/un-hashed identifiers, collisions, and every supported provenance allocation.

Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 91 [x]: Run the exact Topological Torsion environment strict test target added in Step 89.

Step 92 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 93 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 94 [x]: Implement both pinned `getTopologicalTorsionGenerator` overloads as one safe Rust factory over the shared generator core.

Step 95 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 97 [x]: Extend the shared generator JSON dispatch with Topological Torsion arguments, environment generator, and AtomPair invariant generator variants.

Step 98 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 99 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 100 [x]: Extend the shared scalar and bulk fingerprint-type dispatch with Topological Torsion for all four vector forms.

Step 101 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 103 [x]: Add focused strict generator tests for all four scalar outputs, all four ordered bulk outputs, options mutation, JSON roundtrip, custom invariant generator ownership, count bounds, extra bits, chirality, shortest paths, and repeated calls.

Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 105 [x]: Run the exact modern Topological Torsion generator strict test target added in Step 103.

Step 106 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 107 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 108 [x]: Implement pinned legacy `getTopologicalTorsionFingerprint` as an adapter over the shared invariant, path, environment, and torsion-code kernels with the compatibility result size.

Step 109 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 111 [x]: Implement pinned legacy `TorsionFpCalc` as the sole legacy-to-modern hashed-count adapter.

Step 112 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 113 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 114 [x]: Implement pinned legacy `getHashedTopologicalTorsionFingerprint` over `TorsionFpCalc`.

Step 115 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 117 [x]: Implement pinned legacy `getHashedTopologicalTorsionFingerprintAsBitVect` with exact block sizing and four-bit-versus-non-four-bit threshold branches.

Step 118 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 119 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 120 [x]: Add focused strict legacy tests for exact sparse ids/counts, compatibility size, hashed counts, bit thresholds, non-divisible sizes, custom invariants, selections, chirality, invalid parameters, and modern-versus-legacy relationships.

Step 121 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 122 [x]: Run the exact legacy Topological Torsion strict test target added in Step 120.

Step 123 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 125 [x]: Add distinct public Rust Topological Torsion generator, scalar convenience, legacy compatibility, result, parameter, and output-request types without changing the existing RDKFingerprint API.

Step 126 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 127 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 128 [x]: Add focused public Rust API tests for defaults, naming separation, vector forms, provenance access, structured errors, source-molecule immutability, and repeated-call determinism.

Step 129 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 130 [x]: Run the exact public Rust Topological Torsion API strict test target added in Step 128.

Step 131 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 132 [x]: Add the Python Topological Torsion generator, scalar vector forms, legacy compatibility functions, mutable options, shared `AdditionalOutput`, and JSON bindings with pinned wrapper defaults.

Step 133 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 134 [x]: Add the pinned Python atom-code constants, `GetAtomCode`, `ExplainAtomCode`, `pyScorePath`, `ExplainPathScore`, and `GetTopologicalTorsionFingerprintAsIds` helper surface over the shared Rust atom-code and torsion-code core.

Step 135 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 136 [x]: Add focused Python tests for generator defaults, scalar and bulk vector forms, options mutation, empty-list conversion, custom invariants, provenance, JSON, legacy wrappers, helper utilities, typed errors, and Rust-versus-Python determinism.

Step 137 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 138 [x]: Run the exact focused Python Topological Torsion test target added in Step 136.

Step 139 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 140 [x]: Add a pinned deterministic RDKit Topological Torsion oracle generator covering every modern, legacy, provenance, JSON, selection, chirality, path, count, collision, and error profile.

Step 141 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 142 [x]: Generate checksummed focused expected data with source-version metadata and deterministic logs through the pinned oracle generator.

Step 143 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 145 [x]: Add one exact focused Rust parity test consuming every generated modern, legacy, provenance, JSON, and error profile without filtering or tolerance.

Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 147 [x]: Run the exact focused Rust Topological Torsion parity test added in Step 145 with strict operation contracts.

Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 149 [x]: Extend the pinned oracle generator with a source-meaningful Topological Torsion profile matrix over the maintained 5,000-row corpus.

Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 151 [x]: Generate checksummed 5,000-row expected data with deterministic logs through the pinned oracle profile generator.

Step 152 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 153 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 154 [x]: Add an exact 5,000-row Topological Torsion parity test that consumes every generated profile without sampling, fallback, tolerance, or ignored mismatches.

Step 155 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 156 [x]: Run the exact 5,000-row Topological Torsion parity test added in Step 154 in release mode with deterministic parallel execution.

Step 157 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 159 [x]: Add ordered batch Topological Torsion APIs that delegate each record to the single scalar Rust core and preserve indexed structured errors.

Step 160 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 161 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 162 [x]: Add focused strict batch tests for order, scalar equality, invalid records, thread counts, deterministic scheduling, and source-molecule immutability.

Step 163 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 164 [x]: Run the exact Topological Torsion batch strict test target added in Step 162.

Step 165 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 166 [x]: Update the public Topological Torsion documentation surface with Rust and Python examples, four modern outputs, legacy differences, option semantics, provenance, exact-parity boundaries, errors, and the distinction from RDKFingerprint.

Step 167 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 168 [x]: Run `cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen` to generate `python/cosmolkit.pyi` after the Python API changes.

Step 169 [x]: Read `dev/README.md` before editing `crates/cosmolkit-core/` for the next task.

Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 171 [x]: Synchronize core exports, feature inventories, testing contracts, and parity scope by adding a distinct Topological Torsion `FeatureSpec` for the proven boundary.

Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 173 [x]: Audit final source-line closure, single-core architecture, exact parity, errors, determinism, allocation shape, asymptotic complexity, and benchmarks and write `dev/gap_reports/rdkit_topological_torsion_fingerprint_full_port_validation.md`.

Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 175 [x]: Run `cargo fmt --all`.

Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 177 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.

Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 179 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.

Step 180 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 181 [x]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.

Step 182 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 183 [x]: Run `uv sync --group dev` from the repository root.

Step 184 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 185 [x]: Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml`.

Step 186 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 187 [x]: Run `.venv/bin/pytest`.

Step 188 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 189 [x]: Run `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html`.

Step 190 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.

Step 191 [x]: Run `.venv/bin/basedpyright python/tests python/examples`.
