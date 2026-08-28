# RDKit High-Feasibility Descriptor Families Full Source-Port Plan

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

## Source Pin And Scope

- Source and behavior oracle: RDKit `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`, vendored under `third_party/rdkit/`.
- This plan contains only the descriptor families previously assessed as high feasibility: connectivity/Chi, Hall-Kier/Kappa/Phi, Lipinski count extensions, MQN, LabuteASA, SlogP-VSA, and SMR-VSA.
- The Rust API uses project-native `snake_case` names. RDKit names below are source anchors and oracle identifiers, not a second public naming convention.
- Existing exact implementations of HBA/HBD, fraction Csp3, aromatic ring count, rotatable bonds, Crippen atom contributions, TPSA, ring finding, SMARTS/substructure matching, valence, total-H calculation, and stereochemistry are dependencies to reuse and prove, not functions to duplicate.
- The plan excludes EState/EState-VSA, fragment descriptors, BalabanJ, BertzCT, IPC, BCUT, AUTOCORR, all 3D descriptor families, Gasteiger charges, PEOE-VSA, and custom-property VSA.
- Wiener, Zagreb, and Randic are excluded because they are not registered descriptor ports in the pinned RDKit release; a future independent implementation must not be described as an RDKit source port without a pinned source boundary.
- Unsupported input state discovered inside a selected source function must remain a structured unsupported result until its exact source dependency is ported; it must not receive a descriptor-local heuristic.

## Normative Source Call Graphs

```text
Connectivity
  Rust/Python descriptor facade
    -> Chi0 / Chi1 Python-source behavior
    -> calcChi{0..4}{v,n} thin fixed-order functions
      -> calcChiNv / calcChiNn
        -> detail::hkDeltas / detail::nVals
        -> shared findAllPathsOfLengthN atom-path kernel
    -> calcHallKierAlpha
      -> detail::getAlpha
      -> shared PeriodicTable::getRb0
    -> calcKappa{1,2,3} / calcPhi
      -> kappa{1,2,3}Helper
      -> calcHallKierAlpha
      -> shared findAllPathsOfLengthN atom-path kernel

Lipinski and MQN
  Rust/Python descriptor facade
    -> direct atom counts: calcLipinskiHBA/HBD, calcNumHeavyAtoms/Atoms
    -> shared SMARTS matcher: calcNumHBA/HBD/Heteroatoms/AmideBonds
    -> shared SSSR/ring info: ring-class counts, spiro, bridgehead
    -> shared source-backed stereochemistry assignment: stereo-center counts
    -> calcMQNs
      -> shared atom/valence/H/ring state
      -> existing strict rotatable-bond core

Molecular surface
  Rust/Python descriptor facade
    -> calcLabuteASA
      -> getLabuteAtomContribs
        -> shared PeriodicTable::getRb0
        -> descriptor computed-property cache
    -> calcSlogP_VSA / calcSMR_VSA
      -> getLabuteAtomContribs
      -> existing getCrippenAtomContribs
      -> assignContribsToBins
    -> scalar bin projections index the same vector result
```

## Existing Infrastructure Ownership

| Facility | Current owner | Required disposition |
|---|---|---|
| RDKit `findAllPathsOfLengthN` reproduction | Private `find_all_paths_of_length_n()` in `properties/fingerprint.rs` | Extract to one neutral internal subgraph module; fingerprints and connectivity descriptors must call the same implementation. |
| Outer-electron and periodic-table data | `chemistry/valence.rs` and `model::ElementInfo` | Reuse the single table and extend its authoritative record where required. |
| Exact RDKit Rb0 values | Private `add_hs_rdkit_rb0()` in `operations/ops/hydrogens.rs` | Move into the authoritative periodic-table layer; hydrogen placement, Hall-Kier, and LabuteASA must share it. |
| Descriptor computed properties | `model/molecule.rs::ComputedPropertyCache` | Extend this observational, thread-safe, clone-copied, topology-invalidated cache; do not use public molecule properties or operation bypasses. |
| Crippen per-atom contributions | `properties/descriptors.rs` | Reuse one implementation for Crippen and both VSA families; expose only a private contribution view as needed. |
| HBA/HBD and rotatable-bond definitions | `properties/descriptors.rs` | Preserve as the sole existing standard definitions; add distinct Lipinski N/O and N/O-H counts without renaming or replacing them. |
| Ring and ring-family state | `chemistry/rings.rs` and molecule derived cache | Reuse source-backed SSSR ordering and membership; do not add descriptor-local cycle finding. |
| Stereo perception | `chemistry/stereo.rs` | Reuse through a non-mutating clone when source assignment is required; descriptor calls must not mutate the caller. |
| Golden-data framework | `tools/testdata/rdkit/`, `testdata/descriptors/`, and `tools/chembl_parity/` | Extend the existing generator/profile/schema pipeline; do not add a temporary second oracle workflow. |

## Single-Core Target Architecture

- `properties/descriptors.rs` remains the narrow public facade and owns public option/result types and re-exports.
- `properties/descriptors/connectivity.rs` owns connectivity delta generation, Chi, Hall-Kier, Kappa, and Phi only.
- `properties/descriptors/lipinski.rs` owns direct counts, source SMARTS declarations, ring classification, spiro/bridgehead, and stereo-center descriptor projections only.
- `properties/descriptors/mqn.rs` owns one fixed-order `[u32; 42]` MQN implementation; no scalar MQN copies are allowed.
- `properties/descriptors/mol_surface.rs` owns Labute contributions/ASA and one shared VSA binning function; SlogP-VSA and SMR-VSA only select property contributions and bins.
- A neutral chemistry subgraph module owns the sole `findAllPathsOfLengthN` reproduction and its atom-path/ring-path ordering.
- The authoritative periodic-table record owns the sole Rb0 value per element.
- Public fixed-width outputs use `[u32; 42]`, `[f64; 12]`, and `[f64; 10]`; custom-bin VSA calls return `Vec<f64>` because their width is caller-defined.
- Python scalar names for VSA bins are thin index projections over the vector functions, not separately calculated descriptors.
- Read-only descriptor calls do not enter `molecule_ops!`; any logical cache write stays in `ComputedPropertyCache` and follows its clone/invalidation contract.
- Static architecture tests must reject a second path enumerator, Rb0 table, ring classifier, MQN engine, Labute engine, VSA binning loop, or Crippen contribution engine.

## Function-By-Function Source Ledger

### Connectivity

| Source function | Pinned source | Required COSMolKit disposition |
|---|---|---|
| `detail::hkDeltas` | `ConnectivityDescriptors.cpp` | Port atomic-number branches, total-H use, square-root inversion, zero handling, cache lookup/write, and `force`. |
| `detail::nVals` | `ConnectivityDescriptors.cpp` | Port outer-electron minus total-H calculation, square-root inversion, zero handling, cache lookup/write, and `force`. |
| `detail::getAlpha` | `ConnectivityDescriptors.cpp` | Port every atomic-number/hybridization branch and the `found` distinction exactly. |
| `calcChiNv` | `ConnectivityDescriptors.cpp` | Port source path order, multiplication order, and the closed-path rule `if (p[n] != p[0])`. |
| `calcChiNn` | `ConnectivityDescriptors.cpp` | Port the same path and closed-path behavior using nVals. |
| `calcChi0v` | `ConnectivityDescriptors.cpp` | Port ordered accumulation over cached HK deltas. |
| `calcChi1v` | `ConnectivityDescriptors.cpp` | Port source bond iteration and endpoint product order. |
| `calcChi2v` | `ConnectivityDescriptors.cpp` | Thin call to the sole `calcChiNv(order=2)` core. |
| `calcChi3v` | `ConnectivityDescriptors.cpp` | Thin call to the sole `calcChiNv(order=3)` core. |
| `calcChi4v` | `ConnectivityDescriptors.cpp` | Thin call to the sole `calcChiNv(order=4)` core. |
| `calcChi0n` | `ConnectivityDescriptors.cpp` | Port ordered accumulation over cached nVals. |
| `calcChi1n` | `ConnectivityDescriptors.cpp` | Port source bond iteration and endpoint product order. |
| `calcChi2n` | `ConnectivityDescriptors.cpp` | Thin call to the sole `calcChiNn(order=2)` core. |
| `calcChi3n` | `ConnectivityDescriptors.cpp` | Thin call to the sole `calcChiNn(order=3)` core. |
| `calcChi4n` | `ConnectivityDescriptors.cpp` | Thin call to the sole `calcChiNn(order=4)` core. |
| `Chi0` | `rdkit/Chem/GraphDescriptors.py` | Port graph-degree filtering and ordered square-root accumulation independently from valence Chi. |
| `Chi1` | `rdkit/Chem/GraphDescriptors.py` | Port bond endpoint-degree products, zero filtering, and ordered square-root accumulation. |
| `calcHallKierAlpha` | `ConnectivityDescriptors.cpp` | Port dummy skipping, explicit alpha table, Rb0 fallback, optional atom contributions, and source order. |
| `kappa1Helper` | `ConnectivityDescriptors.cpp` | Port denominator and polynomial expression order exactly. |
| `kappa2Helper` | `ConnectivityDescriptors.cpp` | Port denominator and polynomial expression order exactly. |
| `kappa3Helper` | `ConnectivityDescriptors.cpp` | Port odd/even heavy-atom branches and expression order exactly. |
| `calcKappa1` | `ConnectivityDescriptors.cpp` | Port P1/A/alpha inputs and helper call. |
| `calcKappa2` | `ConnectivityDescriptors.cpp` | Port length-2 path count, A/alpha inputs, and helper call. |
| `calcKappa3` | `ConnectivityDescriptors.cpp` | Port length-3 path count, signed A conversion, and helper call. |
| `calcPhi` | `ConnectivityDescriptors.cpp` | Port zero-heavy-atom result and the shared kappa1/kappa2 calculation path. |

### Lipinski Counts And Ring Topology

| Source function | Pinned source | Required COSMolKit disposition |
|---|---|---|
| `calcLipinskiHBA` | `Lipinski.cpp` | Port the direct N/O atom count as a distinct API from the existing SMARTS HBA. |
| `calcLipinskiHBD` | `Lipinski.cpp` | Port the sum of N/O total hydrogens including neighbors exactly. |
| `calcNumHBD` | `Lipinski.cpp` | Reuse the existing SMARTS implementation after a source-and-test proof; do not reimplement. |
| `calcNumHBA` | `Lipinski.cpp` | Reuse the existing recursive-SMARTS implementation after a source-and-test proof; do not reimplement. |
| `calcNumHeteroatoms` | `Lipinski.cpp` | Port the exact `[!#6;!#1]` source SMARTS through the shared matcher. |
| `calcNumAmideBonds` | `Lipinski.cpp` | Port the exact `C(=[O;!R])N` source SMARTS through the shared matcher. |
| `calcNumRings` | `Lipinski.cpp` | Project the source SSSR ring count from shared ring state. |
| `calcNumHeterocycles` | `Lipinski.cpp` | Port atom-ring iteration and first non-carbon short circuit. |
| `calcNumAromaticRings` | `Lipinski.cpp` | Reuse the existing implementation after proving bond-ring and ordering equivalence. |
| `calcNumSaturatedRings` | `Lipinski.cpp` | Port all-single and non-aromatic bond-ring classification. |
| `calcNumAliphaticRings` | `Lipinski.cpp` | Port any-non-aromatic bond-ring classification. |
| `calcNumAromaticHeterocycles` | `Lipinski.cpp` | Port all-aromatic plus any non-carbon endpoint classification. |
| `calcNumAromaticCarbocycles` | `Lipinski.cpp` | Port all-aromatic plus all-carbon endpoint classification. |
| `calcNumAliphaticHeterocycles` | `Lipinski.cpp` | Port any-non-aromatic plus any non-carbon endpoint classification. |
| `calcNumAliphaticCarbocycles` | `Lipinski.cpp` | Port any-non-aromatic plus all-carbon endpoint classification. |
| `calcNumSaturatedHeterocycles` | `Lipinski.cpp` | Port all-single/non-aromatic plus any non-carbon endpoint classification. |
| `calcNumSaturatedCarbocycles` | `Lipinski.cpp` | Port all-single/non-aromatic plus all-carbon endpoint classification. |
| `calcNumSpiroAtoms` | `Lipinski.cpp` | Port pairwise atom-ring intersection, one-atom criterion, uniqueness, and optional ordered atom list. |
| `calcNumBridgeheadAtoms` | `Lipinski.cpp` | Port bond-ring intersection, shared-bond endpoint counting, one-count endpoint selection, uniqueness, and optional ordered atom list. |
| `hasStereoAssigned` | `Lipinski.cpp` | Map source `_StereochemDone` to COSMolKit's authoritative stereo-valid state without inventing a property flag. |
| `numAtomStereoCenters` | `Lipinski.cpp` | Port conditional clone/assignment and `_ChiralityPossible` count without mutating the caller. |
| `numUnspecifiedAtomStereoCenters` | `Lipinski.cpp` | Port conditional clone/assignment and unspecified-tag filter without mutating the caller. |
| `calcNumHeavyAtoms` | `MolDescriptors.cpp` | Port source heavy-atom semantics as a thin graph projection. |
| `calcNumAtoms` | `MolDescriptors.cpp` | Port `getNumAtoms(onlyExplicit=false)` semantics including implicit hydrogens. |
| `calcFractionCSP3` | `Lipinski.cpp` | Reuse the existing exact implementation after a source-and-test proof. |
| `calcNumRotatableBonds` | `Lipinski.cpp` | Reuse the existing option-complete implementation; MQN must call its source-default route. |

### MQN

| Source function | Pinned source | Required COSMolKit disposition |
|---|---|---|
| `calcMQNs` | `MQN.cpp` | Port the complete 42-entry atom, charge, donor/acceptor, degree, bond, aromatic redistribution, rotatable-bond, ring-size, multi-ring atom, and multi-ring bond algorithm in source order. |
| `CalcMQNs` wrapper | `Descriptors/Wrap/rdMolDescriptors.cpp` | Project the fixed `[u32; 42]` core to Python without a second calculation or reordered fields. |

### Labute And VSA

| Source function | Pinned source | Required COSMolKit disposition |
|---|---|---|
| `getLabuteAtomContribs` | `MolSurf.cpp` | Port Rb0 lookup, bond scale factors, aromatic branch, clamped distance, implicit-H aggregate, `1e-4` threshold, expression/summation order, cache lookup/write, and `force`. |
| `calcLabuteASA` | `MolSurf.cpp` | Port cache lookup and delegation while preserving source first-call cache behavior across `includeHs`. |
| `_LabuteHelper` | `rdkit/Chem/MolSurf.py` | Project `[h_contribution, atom_contributions...]` ordering from the same core. |
| `assignContribsToBins` | `MolSurf.cpp` | Port equal-length checks, `upper_bound` boundary behavior, source iteration, and accumulation order once. |
| `calcSlogP_VSA` | `MolSurf.cpp` | Port 11 default boundaries, custom bins, Labute/Crippen force propagation, and 12-bin result through the shared binning core. |
| `calcSMR_VSA` | `MolSurf.cpp` | Port 9 default boundaries, custom bins, Labute/Crippen force propagation, and 10-bin result through the shared binning core. |
| `SlogP_VSA1..12` | `rdkit/Chem/MolSurf.py` | Add thin one-based public projections over the sole 12-bin vector result. |
| `SMR_VSA1..10` | `rdkit/Chem/MolSurf.py` | Add thin one-based public projections over the sole 10-bin vector result. |

## Behavioral Matrices

### Shared Graph And Chemistry State

- Empty, single-atom, disconnected, explicit-H, isotopic-H, dummy-atom, charged, radical, metal, and elements through the complete supported periodic table.
- Unsanitized and sanitized inputs are separate branches; missing source-required valence, ring, hybridization, or stereo state must reproduce the source transition or return the existing structured unsupported/error result.
- Linear, branched, three-membered, aromatic, fused, bridged, spiro, macrocyclic, and multi-ring systems.
- Path lengths 1 through at least 6, open paths, source-recognized closed paths, and the Chi last-element closed-path rule.
- Repeated calls, cloned molecules, mixed descriptor call order, and concurrent reads from the same molecule.

### Cache Semantics

- Cold, warm, and `force=true` calls for HK deltas, nVals, Labute atom/H contributions, LabuteASA, and the existing Crippen contribution cache.
- The source's single Labute cache slot is tested across `include_hs=true/false` call order; `force` is the only recomputation override unless source evidence says otherwise.
- Clone copies current computed values but subsequent writes remain independent; topology edits invalidate descriptor computed properties and coordinate-only edits do not invalidate topology-only descriptor values.
- Cache state is observational and excluded from molecule equality; descriptor calls never mutate public topology, coordinates, properties, or derived chemistry blocks.

### Numeric And Vector Semantics

- Floating-point expression order, path order, bond order, atom order, ring order, and accumulation order match pinned source before bit-level comparison is accepted.
- Every finite floating result is compared by IEEE-754 bits in focused and golden tests; NaN and infinity are compared by explicit classification and sign/payload policy derived from the source branch.
- Integer counts and all 42 MQN entries compare exactly; no aggregate checksum may substitute for per-entry comparison.
- Default and custom VSA bins cover values below, equal to, and above every boundary, repeated boundaries, empty bins, and non-sorted input according to pinned wrapper behavior.
- VSA vectors compare every bin, preserve `upper_bound` equality placement, and verify that bin sums equal the source atom-contribution sum under the same operation order.

## Validation Ladder

1. Focused Rust unit and integration tests reproduce every official pinned RDKit descriptor case plus boundary fixtures absent from upstream tests.
2. `smiles_small` checks every selected function, option, cache branch, atom contribution, and vector entry against committed pinned-RDKit expected data.
3. `smiles_5000` repeats the full source-meaningful matrix without row filtering, mismatch allowlists, similarity scoring, or silent skips.
4. ChEMBL 37 uses the repository's reproducible manifest and streaming parity harness to check every mutually parseable row and every selected profile; summaries record input, eligible, unsupported, compared, and mismatched counts separately.
5. Rust scalar, Rust batch, Python scalar, Python batch/composition, clone, repeated-call, mixed-family, and parallel-read paths all reach the same cores.
6. Benchmarks compare cold/warm cache behavior and representative path/ring sizes against pinned RDKit while source-marker review separately proves asymptotic and allocation equivalence.

## Completion Criteria

- Every source-ledger row has one concrete implementation/reuse proof and focused test artifact, with no approximate, placeholder, heuristic, or silent fallback branch inside the declared boundary.
- Every newly reproduced C/C++ line is copied into the corresponding Rust function with line-level two-axis markers; Python-only source behavior receives equivalent source anchors.
- There is exactly one path enumerator, Rb0 table, ring-state provider, standard HBA/HBD implementation, rotatable-bond implementation, Crippen contribution engine, MQN engine, Labute engine, and VSA binning engine.
- Connectivity cache, Labute cache, clone, invalidation, `force`, and concurrent-read semantics are explicitly verified rather than inferred from scalar output parity.
- All exact counts, atom contribution arrays, MQN entries, VSA bins, and finite floating-point bits match pinned RDKit across focused, small, 5,000-row, and ChEMBL 37 gates.
- The pre-existing strict baseline is recorded before descriptor edits; unrelated baseline failures are neither attributed to this plan nor hidden, and final signoff requires no new failures plus resolution or explicit external blocking disposition for every baseline failure.
- Rust exports, Python exports, generated stubs, examples, public docs, support metadata, the canonical `VALIDATION.md` parity-scope ledger, and `CHANGELOG.md` are updated only after their corresponding evidence exists.
- Final strict core, strict workspace, Python, docs, typing, source-marker, duplicate-implementation, and formatting gates pass.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict` to establish the pre-port strict baseline.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit the baseline result and write every pre-existing failure with ownership and blocking disposition to `dev/gap_reports/rdkit_high_feasibility_descriptors_baseline.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Audit every source-ledger function and dependency against pinned RDKit and current COSMolKit and write exact source anchors, existing owners, gaps, and API dispositions to `dev/gap_reports/rdkit_high_feasibility_descriptors_source_inventory.md`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Modify the descriptor module layout to create the connectivity, Lipinski, MQN, and molecular-surface ownership boundaries defined by this plan without changing behavior.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Add architecture tests that reject duplicate descriptor family cores and preserve every existing descriptor export.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Run the focused descriptor architecture tests.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Extract `find_all_paths_of_length_n` into the neutral shared subgraph module while retaining one source-marked implementation.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Add focused shared-path tests for atom-path ordering, ring paths, disconnected graphs, hydrogens, rooted paths, and existing fingerprint behavior.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Run the focused shared-path and fingerprint regression tests.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Consolidate the RDKit Rb0 table into the authoritative periodic-table implementation and make hydrogen placement consume it.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Add focused periodic-table tests for every Rb0 entry, invalid atomic numbers, and unchanged hydrogen-placement coordinates.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Run the focused Rb0 and hydrogen-placement regression tests.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Extend `ComputedPropertyCache` with typed HK-delta, nVal, and Labute contribution/value entries using the existing clone and invalidation model.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Add focused cache-lifecycle tests for cold, warm, forced, cloned, invalidated, mixed-call-order, and concurrent descriptor reads.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Run the focused descriptor computed-cache lifecycle tests.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port `detail::hkDeltas` completely into the connectivity module.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port `detail::nVals` completely into the connectivity module.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Port `detail::getAlpha` completely into the connectivity module.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Add focused tests for HK deltas, nVals, alpha branches, total-H handling, periodic rows, dummy atoms, caches, and `force`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Run the focused connectivity primitive tests.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Port `calcChiNv` completely through the shared path kernel.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Port `calcChiNn` completely through the shared path kernel.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Port `calcChi0v` completely as the fixed-order valence projection.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Port `calcChi1v` completely as the source bond-order projection.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Port `calcChi2v` completely as a thin call to the generic valence core.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Port `calcChi3v` completely as a thin call to the generic valence core.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Port `calcChi4v` completely as a thin call to the generic valence core.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Port `calcChi0n` completely as the fixed-order nVal projection.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Port `calcChi1n` completely as the source bond-order projection.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Port `calcChi2n` completely as a thin call to the generic nVal core.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Port `calcChi3n` completely as a thin call to the generic nVal core.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Port `calcChi4n` completely as a thin call to the generic nVal core.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Port Python-source `Chi0` completely as the project-native graph-degree descriptor.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Port Python-source `Chi1` completely as the project-native bond-degree descriptor.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Add focused tests for every generic and fixed Chi function including open, closed, three-membered-ring, disconnected, explicit-H, and heavy-element branches.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run the focused Chi descriptor tests.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Port `calcHallKierAlpha` completely with atom-contribution output and shared Rb0 fallback.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Port `kappa1Helper` completely into the private connectivity core.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Port `kappa2Helper` completely into the private connectivity core.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Port `kappa3Helper` completely into the private connectivity core.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Port `calcKappa1` completely through the shared Hall-Kier core.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Port `calcKappa2` completely through the shared path and Hall-Kier cores.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Port `calcKappa3` completely through the shared path and Hall-Kier cores.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Port `calcPhi` completely through the shared helper path.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Add focused tests for every Hall-Kier alpha branch, atom contribution, kappa helper boundary, Kappa1-3 result, and Phi zero-heavy-atom branch.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Run the focused Hall-Kier, Kappa, and Phi tests.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Write source-and-test reuse proofs for existing `calcNumHBA`, `calcNumHBD`, `calcFractionCSP3`, `calcNumAromaticRings`, and `calcNumRotatableBonds` implementations into the descriptor source inventory.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Port `calcLipinskiHBA` completely as the distinct direct N/O count.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Port `calcLipinskiHBD` completely as the distinct N/O total-hydrogen sum.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Port `calcNumHeteroatoms` completely through the shared source SMARTS matcher.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Port `calcNumAmideBonds` completely through the shared source SMARTS matcher.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Port `calcNumHeavyAtoms` completely as a graph projection.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Port `calcNumAtoms` completely with `onlyExplicit=false` semantics.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Add focused tests for direct Lipinski HBA/HBD, standard HBA/HBD reuse, heteroatoms, amide bonds, heavy atoms, and total atoms.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Run the focused Lipinski direct-count and SMARTS-count tests.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Port `calcNumRings` completely through shared SSSR state.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Port `calcNumHeterocycles` completely through shared atom-ring state.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Port `calcNumSaturatedRings` completely through shared bond-ring state.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Port `calcNumAliphaticRings` completely through shared bond-ring state.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Port `calcNumAromaticHeterocycles` completely through shared bond-ring state.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Port `calcNumAromaticCarbocycles` completely through shared bond-ring state.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Port `calcNumAliphaticHeterocycles` completely through shared bond-ring state.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Port `calcNumAliphaticCarbocycles` completely through shared bond-ring state.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Port `calcNumSaturatedHeterocycles` completely through shared bond-ring state.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Port `calcNumSaturatedCarbocycles` completely through shared bond-ring state.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Add focused tests for every ring classification on isolated, aromatic, saturated, partially unsaturated, heterocyclic, carbocyclic, fused, bridged, spiro, and macrocyclic fixtures.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run the focused Lipinski ring-classification tests.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Port `calcNumSpiroAtoms` completely with ordered optional atom output.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Port `calcNumBridgeheadAtoms` completely with ordered optional atom output.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Add focused spiro and bridgehead tests for pairwise ring intersections, duplicate suppression, output ordering, fused systems, and shared-bond endpoints.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Run the focused spiro and bridgehead tests.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Port `hasStereoAssigned` completely onto the authoritative stereo-valid state.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Port `numAtomStereoCenters` completely through source-backed stereochemistry assignment on a conditional clone.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Port `numUnspecifiedAtomStereoCenters` completely through the same conditional-clone path.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Add focused stereo-center count tests for preassigned, unassigned, specified, unspecified, repeated-sanitize, and caller-nonmutation branches.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Run the focused stereo-center descriptor tests.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Port `calcMQNs` completely as the sole fixed-order 42-entry core.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Add focused MQN tests for every index, aromatic redistribution, charge, donor/acceptor, degree, ring-size, multi-ring, and rotatable-bond branch.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Run the focused MQN tests.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Port `getLabuteAtomContribs` completely through shared Rb0 and typed computed-property caches.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Port `calcLabuteASA` completely through the sole contribution core.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Port `_LabuteHelper` completely as the hydrogen-first contribution projection.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Add focused Labute tests for atom/H contributions, all bond types, aromatic bonds, implicit H, cache order, `force`, clone, invalidation, and exact summation bits.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Run the focused Labute contribution and ASA tests.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Port `assignContribsToBins` completely as the sole VSA binning core.
Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Add focused binning tests for exact boundaries, repeated boundaries, empty/custom bins, length errors, non-finite properties, and source accumulation order.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Run the focused shared VSA binning tests.
Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [x]: Port `calcSlogP_VSA` completely through shared Labute, Crippen, and binning cores.
Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [x]: Port `calcSMR_VSA` completely through shared Labute, Crippen, and binning cores.
Step 180 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 181 [x]: Add focused SlogP-VSA and SMR-VSA tests for default/custom bins, every boundary, `force`, cache interactions, all vector bins, and scalar projections.
Step 182 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 183 [x]: Run the focused SlogP-VSA and SMR-VSA tests.
Step 184 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 185 [x]: Add project-native Rust exports for every completed scalar, contribution, fixed-vector, and custom-bin descriptor without exposing implementation helpers.
Step 186 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 187 [x]: Add project-native Python bindings whose scalar and vector projections delegate to the same Rust cores.
Step 188 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 189 [x]: Generate `python/cosmolkit.pyi` with the project stub generator.
Step 190 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 191 [x]: Add Rust and Python API-composition tests covering scalar/vector calls, mixed descriptor order, clones, repeated calls, and parallel reads.
Step 192 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 193 [x]: Run the focused Rust and Python descriptor API-composition tests.
Step 194 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 195 [x]: Extend the unified pinned-RDKit descriptor golden generator and schema with every selected scalar, contribution array, cache profile, MQN entry, and VSA bin.
Step 196 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 197 [x]: Run the focused tests for the descriptor golden generator and schema.
Step 198 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 199 [x]: Add focused parity fixtures covering every behavioral-matrix branch and every official pinned RDKit descriptor regression used by the selected families.
Step 200 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 201 [x]: Run the focused pinned-RDKit descriptor parity tests.
Step 202 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 203 [x]: Generate the complete `smiles_small` expected-data profile through `tools/testdata/rdkit/generate_all.py`.
Step 204 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 205 [x]: Run the complete `smiles_small` high-feasibility descriptor parity suite in strict release mode.
Step 206 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 207 [x]: Generate the complete `smiles_5000` expected-data profile through `tools/testdata/rdkit/generate_all.py`.
Step 208 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 209 [x]: Run the complete `smiles_5000` high-feasibility descriptor parity suite in strict release mode.
Step 210 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 211 [x]: Extend the reproducible ChEMBL 37 streaming harness with the complete selected descriptor matrix and explicit record/observation accounting.
Step 212 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 213 [x]: Run the complete ChEMBL 37 high-feasibility descriptor parity phase against pinned RDKit.
Step 214 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 215 [x]: Write the focused, small, 5,000-row, and ChEMBL 37 record/observation/mismatch results to `dev/gap_reports/rdkit_high_feasibility_descriptors_full_port_validation.md` without historical rerun language.
Step 216 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 217 [x]: Add cold-cache, warm-cache, path-length, ring-complexity, MQN, Labute, and VSA benchmarks against comparable pinned-RDKit calls.
Step 218 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 219 [x]: Run the high-feasibility descriptor benchmark suite and record unexplained performance or complexity gaps in the validation report.
Step 220 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 221 [x]: Audit every selected source function for line-level two-axis markers and write any unresolved behavioral or complexity marker to the validation report.
Step 222 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 223 [x]: Audit the final architecture for duplicate algorithms, descriptor-local chemistry, cache bypasses, mutable public storage, and alternative API paths and write the result to the validation report.
Step 224 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 225 [x]: Correct and strengthen the descriptor architecture guard so it checks the actual VSA binning core and the selected public facade delegation paths.
Step 226 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 227 [x]: Run the focused descriptor architecture test.
Step 228 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 229 [x]: Port one authoritative source-backed ring-state acquisition path and route every selected Lipinski ring projection, reused aromatic-ring projection, spiro/bridgehead function, and MQN through it without changing source ring ordering.
Step 230 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 231 [x]: Add focused ring-state tests covering initialized reuse, cold initialization, clone behavior, topology invalidation, descriptor call order, spiro/bridgehead, aromatic-ring reuse, and MQN.
Step 232 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 233 [x]: Run the focused ring-state descriptor tests.
Step 234 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 235 [x]: Port RDKit `CrippenParamCollection` parsed-query reuse and typed per-atom Crippen contribution caching into the single existing Crippen core with exact `force`, clone, and topology-invalidation behavior.
Step 236 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 237 [x]: Add focused Crippen/VSA tests for parsed-query singleton reuse, cold/warm/forced contribution calls, clone, topology invalidation, mixed LogP/MR/VSA call order, and parallel reads.
Step 238 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 239 [x]: Run the focused Crippen and molecular-surface cache tests.
Step 240 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 241 [x]: Audit and source-align the remaining MQN hot path, including total-H access, ring membership, HBA/HBD reuse, and rotatable-bond dispatch, without adding descriptor-local chemistry or heuristics.
Step 242 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 243 [x]: Run the focused MQN tests and the maintained 5,000-row high-feasibility descriptor parity test.
Step 244 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 245 [x]: Rerun the isolated descriptor benchmark and record the post-cache ring, MQN, Crippen/VSA, Labute, and Chi ratios.
Step 246 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 247 [x]: Audit and port source-equivalent corrections for any remaining material Chi or Labute call-path, allocation, ownership, cache-copy, or wrapper cost found by the repeated benchmark.
Step 248 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 249 [x]: Run focused Chi, Labute, VSA, ring, and MQN tests plus the maintained 5,000-row high-feasibility descriptor parity test.
Step 250 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 251 [x]: Rerun the final isolated descriptor benchmark, update every affected two-axis marker, and record the final performance disposition in the validation report.
Step 252 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 253 [x]: Update Rust docs, Python docs, examples, support metadata, and the descriptor inventory to describe only the validated completed surface.
Step 254 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 255 [x]: Update the canonical `VALIDATION.md` parity-scope ledger with the final complete evidence table and comparison details; do not recreate the removed duplicate `dev/parity_scope.md` ledger.
Step 256 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 257 [x]: Update the current unreleased `CHANGELOG.md` section with the completed descriptor capabilities and verified corpus scope.
Step 258 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 259 [x]: Run `cargo fmt --all`.
Step 260 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 261 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 262 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 263 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 264 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 265 [x]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.
Step 266 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 267 [x]: Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml`.
Step 268 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 269 [x]: Run `.venv/bin/pytest`.
Step 270 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 271 [x]: Run `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html`.
Step 272 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 273 [x]: Run `.venv/bin/basedpyright python/tests python/examples`.
Step 274 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 275 [x]: Audit the final baseline, validation, marker, architecture, API, docs, and test artifacts against every completion criterion and write the signoff decision to `dev/gap_reports/rdkit_high_feasibility_descriptors_full_port_validation.md`.
