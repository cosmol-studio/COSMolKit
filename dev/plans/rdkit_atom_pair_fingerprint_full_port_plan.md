# RDKit AtomPair Fingerprint Full Source-Port Plan

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

## Source Pin And Selected Boundary

- Chemistry oracle and source tree: RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`, vendored under `third_party/rdkit/`.
- The normative implementation is the modern `FingerprintGenerator` path rooted at `RDKit::AtomPair::getAtomPairGenerator()`; deprecated `AtomPairs.cpp`, `rdMolDescriptors.cpp`, and `rdkit/Chem/AtomPairs/Pairs.py` entry points are compatibility behavior only.
- The selected output boundary includes exact sparse-count, folded count, sparse-bit, and explicit-bit fingerprints; count simulation; custom atom invariants; chirality; 2D and 3D distance modes; conformer selection; `fromAtoms`; `ignoreAtoms`; every AtomPair-supported `AdditionalOutput`; generator metadata/JSON; and ordered parallel batch generation.
- The project-facing API remains Rust-native and method-oriented; RDKit-style names are source anchors and oracle labels, not a reason to create a second public naming system.
- Topological-torsion fingerprint generation is outside this plan, but `AtomPairAtomInvGenerator`'s source-defined topological-torsion correction remains part of the shared invariant-generator port so a future torsion port will not need a duplicate atom-code implementation.
- `AtomPairs.cpp::setAtomPairBit()` and its `updateElement()` helpers are unreachable from the pinned modern AtomPair path and from the current deprecated facade, so they are source-history evidence rather than a second implementation target.
- `Pairs.py::pyScorePair()`, `ExplainPairScore()`, `Utils.py::ExplainAtomCode()`, and `Utils.py::BitsInCommon()` are diagnostic legacy utilities outside the fingerprint-generation call graph; their encoding expectations are covered at the lower-level code tests without creating parallel production algorithms.

## Normative Call Graph

```text
Rust/Python project API
  -> AtomPairFingerprintGenerator factory
     -> AtomPairArguments
     -> AtomPairAtomInvGenerator
     -> AtomPairEnvGenerator
     -> shared FingerprintGenerator core
        -> stereochemistry/CIP preparation when requested
        -> AtomPairAtomInvGenerator::getAtomInvariants
           -> AtomPairs::getAtomCode
              -> numPiElectrons
              -> existing CIPLabeler::assignCIPLabels when required
        -> AtomPairEnvGenerator::getEnvironments
           -> shared MolOps::getDistanceMat reproduction for use2D=true
           -> shared MolOps::get3DDistanceMat reproduction for use2D=false
           -> AtomPairAtomEnv construction in source pair order
        -> AtomPairAtomEnv::getBitId
           -> AtomPairs::getAtomPairCode for unfolded sparse-count output
           -> existing gboost hash_combine for folded/hashed outputs
        -> existing single result projector
           -> SparseCountFingerprint
           -> folded SparseCountFingerprint
           -> SparseBitFingerprint with count simulation
           -> Fingerprint with count simulation
           -> AdditionalOutput duplication after count simulation
  -> project-native scalar and batch facades delegate to the same generator
  -> deprecated-RDKit semantic adapters, where retained for testing, delegate to the same generator
```

## Existing COSMolKit Infrastructure Ownership

| Existing facility | Current location | Plan disposition |
|---|---|---|
| `AdditionalOutput` allocation/reset model | `crates/cosmolkit-core/src/properties/fingerprint.rs` | Reuse unchanged as the only provenance container; AtomPair only supplies its family-specific update function. |
| `FingerprintArguments` and `FingerprintFuncArguments` | `crates/cosmolkit-core/src/properties/fingerprint.rs` | Reuse as the only common option and call-argument model. |
| `Fingerprint`, `SparseBitFingerprint`, `SparseCountFingerprint` | `crates/cosmolkit-core/src/properties/fingerprint.rs` | Reuse as the only bit/count vector representations; do not add AtomPair-specific vectors. |
| `hash_combine()` | `crates/cosmolkit-core/src/properties/fingerprint.rs` | Reuse as the single 32-bit gboost hash-combine implementation. |
| `RdkitFingerprintMtRng` | `crates/cosmolkit-core/src/properties/fingerprint.rs` | Reuse for generic extra-bits-per-feature behavior; AtomPair defaults to one bit per feature. |
| `duplicate_additional_output_bit()` and `setup_temp_additional_output()` | `crates/cosmolkit-core/src/properties/fingerprint.rs` | Reuse as the only count-simulation provenance projector. |
| `getFingerprintHelper()`, four output projections, and batch projection | Currently methods on `MorganFingerprintGenerator` | Extract once into a family-neutral internal generator core, then make Morgan and AtomPair delegate to it. |
| Dead `build_fingerprint()`/`compute_initial_invariants()` path | `crates/cosmolkit-core/src/properties/fingerprint.rs` | Remove after proving it has no production caller; it must not become an AtomPair shortcut. |
| RDKit `getDistanceMat()` reproduction | Private in `crates/cosmolkit-core/src/chemistry/distgeom.rs` | Move to one shared internal matrix module and keep distgeom on the same implementation. |
| 3D conformer storage | `Molecule::conformers_3d()` and `Conformer3D::coordinates()` | Reuse; add one shared `get3DDistanceMat()` reproduction rather than an AtomPair-local matrix. |
| CIP labeler | `crates/cosmolkit-core/src/chemistry/ciplabeler.rs` | Reuse; do not add an AtomPair-specific CIP algorithm. |
| Batch scheduling/order preservation | `crates/cosmolkit-core/src/properties/batch.rs` | Reuse; add AtomPair callbacks only. |
| Golden generation and corpus manifests | `tools/testdata/rdkit/`, `testdata/`, and `tools/chembl_parity/` | Extend the existing schema and entry points; do not add a standalone temporary generator workflow. |

## Single-Core Target Architecture

- `properties/fingerprint/generator.rs` owns the only common environment-to-vector pipeline, common argument validation, count simulation, random extra-bit expansion, `AdditionalOutput` reset/duplication, JSON envelope, and ordered batch projection.
- `properties/fingerprint/atom_pair.rs` owns only AtomPair constants, atom invariant generation, pair environments, pair code/hash selection, AtomPair arguments, and thin factories/facades.
- Morgan continues to own Morgan invariant/environment behavior but delegates all result projection to `generator.rs`; no output-family algorithm remains duplicated in Morgan helpers.
- A shared chemistry matrix module owns both source-backed 2D topological and 3D Euclidean distance matrices; distgeom and AtomPair consume it without local copies.
- Every scalar Rust method, Python method, batch method, legacy-semantic test adapter, and JSON restoration route constructs or restores the same `AtomPairFingerprintGenerator` and reaches the same `get_fingerprint_helper()` implementation.
- Static architecture tests reject a second AtomPair pair-enumeration loop, a second count-simulation projector, a second Boost hash-combine implementation, and a second fingerprint vector type.

## Function-By-Function Source Ledger

| Source function | Pinned source | Required COSMolKit action |
|---|---|---|
| `RDKit::numPiElectrons` | `Code/GraphMol/Atom.cpp:934` | Port once as a shared atom-property helper using existing valence, hybridization, explicit-H, bond-valence, and aromatic state. |
| `AtomPairs::getAtomCode` | `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:45` | Port atomic-type bucketing, degree subtraction, modulo semantics, pi-electron encoding, conditional CIP assignment, R/S bits, and width postcondition. |
| `AtomPairs::getAtomPairCode` | `Code/GraphMol/Fingerprints/FingerprintUtil.cpp:99` | Port distance precondition, endpoint-order canonicalization, chirality-dependent shift width, and exact `u32` packing. |
| `AtomPairAtomInvGenerator::AtomPairAtomInvGenerator` | `Code/GraphMol/Fingerprints/AtomPairGenerator.cpp:26` | Port both stored flags without inventing family-specific defaults. |
| `AtomPairAtomInvGenerator::getAtomInvariants` | `AtomPairGenerator.cpp:31` | Port source atom iteration/index placement and wrapping-width behavior for the torsion correction branch. |
| `AtomPairAtomInvGenerator::infoString` | `AtomPairGenerator.cpp:45` | Port the exact source metadata string through the common metadata interface. |
| `AtomPairAtomInvGenerator::toJSON` | `AtomPairGenerator.cpp:50` | Port type and both flag fields into the common JSON envelope. |
| `AtomPairAtomInvGenerator::fromJSON` | `AtomPairGenerator.cpp:56` | Port source default-preserving field restoration. |
| `AtomPairAtomInvGenerator::clone` | `AtomPairGenerator.cpp:63` | Port an owned logical clone without sharing mutable generator state. |
| `AtomPairEnvGenerator::getResultSize` | `AtomPairGenerator.cpp:69` | Port exact unfolded result width for chiral and non-chiral codes with checked Rust shift representation. |
| `AtomPairArguments::AtomPairArguments` | `AtomPairGenerator.cpp:77` | Port common defaults, `{1,2,4,8}` bounds, one bit per feature, 2048-bit default, 2D default, and `minDistance <= maxDistance`. |
| `AtomPairArguments::infoString` | `AtomPairGenerator.cpp:89` | Port the exact family-specific metadata string. |
| `AtomPairArguments::toJSON` | `AtomPairGenerator.cpp:94` | Port AtomPair fields followed by the existing common argument fields. |
| `AtomPairArguments::fromJSON` | `AtomPairGenerator.cpp:101` | Port source default-preserving AtomPair and common-field restoration. |
| `AtomPairAtomEnv::updateAdditionalOutput` | `AtomPairGenerator.cpp:109` | Port `bitInfoMap`, duplicate-preserving `atomToBits`, per-endpoint `atomCounts`, and ordered two-atom `atomsPerBit`. |
| `AtomPairAtomEnv::getBitId` | `AtomPairGenerator.cpp:131` | Port custom-invariant modulo by `codeSizeLimit`, canonical endpoint ordering, three-stage gboost hashing, and unfolded code selection. |
| `AtomPairAtomEnv::AtomPairAtomEnv` | `AtomPairGenerator.cpp:170` | Port the two endpoint identifiers and floored distance as the complete immutable environment state. |
| `AtomPairEnvGenerator::getEnvironments` | `AtomPairGenerator.cpp:178` | Port 2D/3D matrix selection, conformer selection, `i < j` order, ignore precedence, rooted-pair inclusion, `floor()` conversion, disconnected behavior, and inclusive distance bounds. |
| `AtomPairEnvGenerator::infoString` | `AtomPairGenerator.cpp:240` | Port the exact environment-generator metadata string. |
| `AtomPairEnvGenerator::toJSON` | `AtomPairGenerator.cpp:245` | Port the environment type discriminator through the common JSON envelope. |
| `getAtomPairGenerator(args, invgen, owns)` | `AtomPairGenerator.cpp:252` | Port default invariant-generator construction and ownership semantics into Rust ownership without a parallel engine. |
| `getAtomPairGenerator(parameters...)` | `AtomPairGenerator.cpp:271` | Port parameter forwarding and default construction as a thin call to the preceding factory. |
| `FingerprintGenerator::getFingerprintHelper` | `FingerprintGenerator.cpp:325` | Reuse the existing port after extracting it from Morgan ownership; close the currently unsupported missing-stereochemistry branch for all generator families. |
| `duplicateAdditionalOutputBit` | `FingerprintGenerator.cpp:438` | Reuse the existing single implementation and prove AtomPair count-simulation provenance passes through it. |
| `setupTempAdditionalOutput` | `FingerprintGenerator.cpp:481` | Reuse the existing single implementation and prove every AtomPair-supported output allocation is preserved. |
| `getSparseCountFingerprint` | `FingerprintGenerator.cpp:505` | Reuse the common projector for exact unfolded pair-code counts. |
| `getSparseFingerprint` | `FingerprintGenerator.cpp:517` | Reuse the common projector for sparse bits and count simulation with the source `u32::MAX` size cap. |
| `getCountFingerprint` | `FingerprintGenerator.cpp:577` | Reuse the common projector for gboost-hashed folded counts. |
| `getFingerprint` | `FingerprintGenerator.cpp:593` | Reuse the common projector for explicit bits, bounds validation, count simulation, and provenance duplication. |
| `mtgetFingerprints` and four batch methods | `FingerprintGenerator.cpp:656-754` | Reuse COSMolKit batch scheduling while preserving source scalar equivalence, order, determinism, and error indexing. |
| `FingerprintGenerator::{infoString,toJSON,fromJSON}` and `generatorFromJSON` | `FingerprintGenerator.cpp:162-317` | Extend the existing common metadata/JSON route with AtomPair type discriminators; do not create AtomPair-only serialization. |
| `getAtomPairFingerprintInternal` | `AtomPairs.cpp:51` | Represent only as a thin compatibility-semantic adapter to the modern generator for exact-vs-folded count selection. |
| `getAtomPairFingerprint` overloads | `AtomPairs.cpp:81,91` | Represent only as thin default/argument adapters to the same sparse-count core. |
| `getHashedAtomPairFingerprint` | `AtomPairs.cpp:105` | Represent only as a thin folded-count adapter to the same core. |
| `getHashedAtomPairFingerprintAsBitVect` | `AtomPairs.cpp:122` | Represent its `nBitsPerEntry == 4` bounds and non-4 linear thresholds using the common count projection without a second pair loop. |
| `AtomPairWrapper::getAtomPairGenerator` | `Wrap/AtomPairWrapper.cpp:23` | Map project-native Python keywords to the same Rust factory and copy user count bounds/invariant configuration exactly once. |
| `AtomPairWrapper::getAtomPairAtomInvGen` | `Wrap/AtomPairWrapper.cpp:45` | Expose the project-native invariant selection through the same Rust invariant-generator type. |
| `rdMolDescriptors` AtomPair wrappers | `Descriptors/Wrap/rdMolDescriptors.cpp:213-330,969-1039` | Use as argument-validation and default-value oracle; do not reproduce RDKit naming as a second Python API. |
| `Pairs.py::GetAtomPairFingerprintAsBitVect` | `rdkit/Chem/AtomPairs/Pairs.py:123` | Cover its unfolded-presence semantics through the common sparse-bit output and do not retain a Python-level conversion loop. |

## Behavioral Matrix Required For Completion

- Atom encoding: every source atomic-number bucket including the fallback bucket, aromatic/non-aromatic state, SP3/non-SP3 state, dative/zero-valence bonds, explicit hydrogens, branch subtraction, degree modulo 7, pi count modulo 3, R/S/no-CIP, and legacy/new stereo-perception state.
- Pair encoding: endpoint reversal, distance 0/1/30/31 boundaries, chirality-dependent shifts, `u32` width, custom invariant values around 511/2047 and `u32::MAX`, and exact-vs-hashed output differences.
- Pair enumeration: empty/single-atom molecules, chains, branches, rings, fused rings, disconnected fragments, explicit H, `fromAtoms`, `ignoreAtoms`, both lists together, duplicates, invalid indices, min/max equality, and inclusive endpoints.
- 3D mode: default and explicit conformer IDs, multiple conformers, exact integer distances, fractional values on both sides of an integer, missing conformers, coordinate-count mismatch, and non-finite coordinates according to the shared coordinate contract.
- Output families: sparse count, folded count, sparse bit, explicit bit, default and custom count bounds, count-simulation on/off, collisions, every allocated `AdditionalOutput`, and source ordering after duplication.
- API combinations: scalar/batch, Rust/Python, JSON round-trip, custom invariants with chirality, 2D/3D, rooted/ignored pairs, output types, repeated calls, mixed fingerprint-family call order, and concurrent calls on shared molecule values.
- Corpora: focused project fixtures, every active row of `smiles_small`, every active row of `smiles_5000`, and every valid ChEMBL 37 row for source-meaningful profiles; no row filtering, mismatch allowlist, similarity-only comparison, or silent skip is permitted.

## Completion Criteria

- Every function in the source ledger has a disposition backed by a concrete code/test artifact or an explicit unreachable/out-of-scope proof.
- Every reproduced C++ line in AtomPair-specific and newly shared dependency functions carries a line-level two-axis marker; reused generic functions keep one authoritative source frame.
- There is exactly one implementation of pair enumeration, pair encoding, gboost hash combine, generic result projection, count simulation, `AdditionalOutput` duplication, 2D distance matrices, 3D distance matrices, and each vector representation.
- Morgan's existing exact-bit and `AdditionalOutput` behavior remains unchanged after generic-core extraction.
- All four modern output forms, all selected options, metadata/JSON, scalar/batch routes, and compatibility-semantic adapters match pinned RDKit exactly on the supported state space.
- Focused, small-corpus, 5,000-row, and ChEMBL 37 audits report zero unexplained mismatches with reproducible manifests and exact comparison counts.
- Benchmarks show no unexplained asymptotic, allocation, or hot-path regression relative to the pinned generator for comparable inputs.
- Rust exports, Python methods, stubs, docs, examples, support metadata, parity scope, and changelog describe the same completed capability only after the evidence exists.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit every function and dependency in the source ledger against the pinned RDKit tree and write `dev/gap_reports/rdkit_atom_pair_fingerprint_source_inventory.md` with exact source ranges, reachability, current COSMolKit ownership, planned target ownership, and a line-coverage ledger.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Modify the existing Morgan-bound generic generator methods into one family-neutral internal `FingerprintGenerator` core and reconnect Morgan to that core without changing public behavior.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Remove the unreachable alternative `build_fingerprint()` projection path and its exclusively owned helpers after recording their caller proof in the source inventory.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Add focused architecture and Morgan regression tests proving that all output forms use the shared projector and that extraction preserves exact bits, counts, count simulation, JSON, and every `AdditionalOutput` field.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::tests`.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Modify the existing RDKit `MolOps::getDistanceMat()` reproduction into one shared internal chemistry-matrix function and reconnect distgeom to it.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Add focused shared-matrix tests for source ordering, disconnected components, empty molecules, weighted-option rejection at the selected boundary, and unchanged distgeom consumers.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::matrices`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port `MolOps::get3DDistanceMat()` from `Code/GraphMol/Matrices.cpp:392` into the shared chemistry-matrix module with exact conformer selection, coordinate order, Euclidean arithmetic, and source error boundaries.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Add focused shared 3D-matrix tests for default and explicit conformers, multiple conformers, empty molecules, fractional coordinates, missing conformers, and malformed coordinate state.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chemistry::matrices::tests::rdkit_3d_distance`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Port `numPiElectrons()` from `Code/GraphMol/Atom.cpp:934` as one shared atom-property helper using the existing source-backed valence and bond-contribution machinery.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Add focused exact-value tests for every `numPiElectrons()` branch including aromatic atoms, SP3 atoms, multiple bonds, dative bonds, zero bonds, explicit hydrogens, and invalid explicit-valence state.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict source_port__atom__num_pi_electrons__line_934`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Implement the AtomPair constant table and integer-width definitions from `FingerprintUtil.h:27-41` in the sole AtomPair module.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port `getAtomCode()` from `FingerprintUtil.cpp:45` with exact branch, pi, element-bucket, chirality, CIP, modulo, and width behavior.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Add focused exact-code tests for every AtomPair atom-code branch and boundary in the behavioral matrix.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::atom_code`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Port `getAtomPairCode()` from `FingerprintUtil.cpp:99` with exact endpoint canonicalization, distance precondition, chirality shift, and `u32` packing.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Add focused exact pair-code tests for reversed endpoints, chirality modes, distance boundaries, maximum codes, and source precondition failures.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::pair_code`.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Port `AtomPairArguments::{AtomPairArguments,infoString,toJSON,fromJSON}` as one complete source-owned argument state machine layered on the existing `FingerprintArguments`.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Add focused argument tests for all defaults, min/max validation, default-preserving partial JSON, full JSON round-trip, invalid JSON types, and exact metadata strings.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::arguments`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Port `AtomPairAtomInvGenerator::{AtomPairAtomInvGenerator,getAtomInvariants,infoString,toJSON,fromJSON,clone}` as one complete invariant-generator state machine using only the shared atom-code helper.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Add focused invariant-generator tests for atom-index ordering, chirality, topological-torsion correction, unsigned boundaries, clone independence, exact metadata, and JSON round-trip.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::atom_invariant_generator`.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Port `AtomPairAtomEnv::AtomPairAtomEnv()` as the sole immutable pair-environment constructor.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Port `AtomPairAtomEnv::getBitId()` with exact custom-invariant modulo, endpoint ordering, shared gboost hashing, and unfolded code dispatch.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Add focused raw-bit-ID tests for exact and hashed modes, chiral widths, collisions, custom invariant boundaries, invalid invariant lengths, and endpoint reversal.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::environment_bit_id`.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Port `AtomPairAtomEnv::updateAdditionalOutput()` into the existing `AdditionalOutput` container with exact duplicate and ordering behavior.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Add focused tests for every isolated and combined AtomPair `AdditionalOutput` allocation before and after hash collisions.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::additional_output`.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Port `AtomPairEnvGenerator::getEnvironments()` using only the shared 2D/3D matrix functions and the sole AtomPair environment constructor.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Add focused ordered-environment tests for every pair-enumeration, distance, conformer, rooted, ignored, disconnected, and bounds branch in the behavioral matrix.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::environment_generation`.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Port `AtomPairEnvGenerator::{getResultSize,infoString,toJSON}` as one complete source-owned environment-generator metadata unit.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Add focused result-size and environment-metadata tests for chiral/non-chiral widths and JSON restoration dispatch.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::environment_generator_metadata`.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Port both `getAtomPairGenerator()` overloads as thin constructors of the shared generic core with Rust-owned invariant-generator semantics.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Extend the existing generic `infoString`, `toJSON`, `fromJSON`, and `generatorFromJSON` routes with AtomPair discriminators that restore the same generator type.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Add focused factory and generator-serialization tests for defaults, custom invariants, ownership independence, exact info strings, partial JSON, full JSON, and all four restored output forms.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::atom_pair::tests::generator_factory`.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Modify the shared fingerprint stereochemistry preparation branch to reproduce RDKit's clone-and-assign behavior by delegating to existing source-backed stereo and CIP components without mutating the input molecule.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Add focused cross-family chirality tests for absent/present `_StereochemDone`, absent/present `_CIPComputed`, R/S/no-label atoms, custom invariants, source immutability, and Morgan non-regression.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::fingerprint::tests::generator_chirality_preparation`.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Implement the project-native Rust AtomPair parameter, result, generator, molecule-method, and crate-export surface as thin owners or delegates of the sole core implementation.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Add public Rust API tests for defaults, all output forms, typed errors, value semantics, repeated-call determinism, mixed fingerprint-family ordering, and absence of RDKit-style duplicate APIs.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test atom_pair_fingerprint_public_api`.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Add exact integration tests for AtomPair sparse-count, folded-count, sparse-bit, and explicit-bit outputs across count-simulation, custom-bound, collision, provenance, and option-interaction branches.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_atom_pair_fingerprint_parity`.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Implement the deprecated `AtomPairs.cpp` sparse-count, folded-count, and count-simulated-bit semantics as private testable adapters that delegate to the same modern generator core.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Add exact compatibility tests proving every retained adapter equals its modern-core equivalent for defaults, `nBitsPerEntry` equal/not-equal to four, custom invariants, chirality, rooted/ignored pairs, and 2D/3D modes.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_atom_pair_compatibility_parity`.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Implement ordered parallel AtomPair batch methods by delegating each record to the scalar Rust API through the existing batch runtime.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Add batch tests for scalar equivalence, record order, invalid-record indexing, thread counts, determinism, shared-input safety, and all public AtomPair result forms selected for batch exposure.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test atom_pair_fingerprint_batch`.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Implement project-native scalar Python AtomPair methods and typed result wrappers that delegate directly to the Rust public API.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Build the Python extension with `.venv/bin/maturin develop --release --manifest-path python/Cargo.toml`.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Add focused scalar Python tests for keyword defaults, result conversion, exact bits/counts/provenance, typed exceptions, immutability, and no duplicate RDKit-style facade names.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Run `.venv/bin/python -m pytest python/tests/test_atom_pair_fingerprint.py -q`.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Implement project-native Python AtomPair batch methods by delegating to the Rust batch surface.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Build the Python extension with `.venv/bin/maturin develop --uv --release --manifest-path python/Cargo.toml` in the repository `.venv`.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Add focused Python batch tests for scalar equivalence, stable order, `n_jobs`, invalid records, progress configuration, repeated calls, and concurrent mixed-family execution.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Run `.venv/bin/python -m pytest python/tests/test_atom_pair_fingerprint_batch.py -q`.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Add the AtomPair expected-data schema, focused profile, corpus profiles, source pin, manifest identity fields, and exact comparison fields to the existing RDKit golden framework.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Implement the deterministic AtomPair oracle generator inside `tools/testdata/rdkit/` with the existing `generate_all.py` entry point and no temporary standalone workflow.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Add generator-contract tests for schema completeness, option identity, source version, record counts, checksums, deterministic ordering, and rejection of incomplete cached output.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Run `.venv/bin/python -m pytest tools/testdata/rdkit/test_expected_schema.py tools/testdata/rdkit/test_atom_pair_generator.py -q`.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Generate the focused and `smiles_small` AtomPair expected-data profiles with the pinned RDKit environment through `tools/testdata/rdkit/generate_all.py`.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Add a public-boundary AtomPair parity test that consumes every focused and `smiles_small` record and compares all configured intermediate and final fields exactly.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_atom_pair_fingerprint_golden -- focused_and_small`.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Generate the complete `smiles_5000` AtomPair expected-data profile with the pinned RDKit environment through `tools/testdata/rdkit/generate_all.py`.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_atom_pair_fingerprint_golden -- smiles_5000`.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Add AtomPair source-meaningful profiles and exact intermediate/final comparison fields to the repository-owned ChEMBL 37 parity workflow.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Add ChEMBL workflow tests for AtomPair shard identity, phase completeness, comparison counts, restart safety, and rejection of filtered or partial results.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Run `.venv/bin/python -m pytest dev/tools/chembl_parity/tests -q`.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Run the complete ChEMBL 37 AtomPair parity phase through the repository-owned workflow with every prepared shard and profile.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Run deterministic AtomPair scalar and parallel benchmarks against pinned RDKit for small, medium, large, 2D, 3D, sparse-count, folded-count, sparse-bit, and explicit-bit profiles.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Audit source-line closure, single-core ownership, exact parity, structured errors, mutation behavior, determinism, allocation shape, asymptotic complexity, benchmark results, and corpus completeness and write `dev/gap_reports/rdkit_atom_pair_fingerprint_full_port_validation.md` with links to ignored machine artifacts.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Port the selected RDKit `getDistanceMat` and `get3DDistanceMat` logically read-only cache behavior into the shared molecule computed-property cache with source-exact 2D and resolved-conformer keys and operation-driven topology and coordinate invalidation.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Add focused tests for cold and warm distance-matrix reuse, clone isolation, topology and coordinate invalidation, resolved conformer keying, parallel initialization, and unchanged exact AtomPair outputs.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Run the focused shared-matrix and AtomPair cache tests under release `op-contracts-strict`.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Rebuild the Python extension in release mode for post-cache benchmarking.
Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Rerun the deterministic AtomPair scalar and parallel benchmark against pinned RDKit and retain its ignored machine artifact.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Update `dev/gap_reports/rdkit_atom_pair_fingerprint_full_port_validation.md` with post-port source closure, targeted results, benchmark evidence, and the final unexplained-gap verdict.
Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [x]: Update support metadata to `supported_with_rdkit_parity` only if the validation report records zero unexplained mismatches and every completion criterion is satisfied.
Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [x]: Update current public Rust docs, Python docs, stubs, examples, README feature summaries, and `VALIDATION.md` to describe the same project-native AtomPair surface and verified corpus evidence without editing historical plans.
Step 180 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 181 [x]: Update the pending release section of `CHANGELOG.md` once with the completed AtomPair implementation and final exact comparison totals.
Step 182 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 183 [x]: Run `cargo fmt --all -- --check`.
Step 184 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 185 [x]: Run `cargo clippy --workspace --all-targets --features cosmolkit-core/op-contracts-strict`.
Step 186 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 187 [ ]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.
Step 188 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 189 [ ]: Run the complete Python test suite against the release-built extension.
Step 190 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 191 [ ]: Run the Sphinx documentation build and doctest suite with warnings treated as errors.
Step 192 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 193 [ ]: Run repository source-marker, duplicate-source-port, generated-stub, support-matrix, expected-data-manifest, and temporary-artifact guards.
Step 194 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 195 [ ]: Move this fully completed plan to `dev/archive/plans/` without rewriting its checked execution history.
Step 196 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 197 [ ]: Update the active and archived plan indexes to reflect the completed plan's final location.
