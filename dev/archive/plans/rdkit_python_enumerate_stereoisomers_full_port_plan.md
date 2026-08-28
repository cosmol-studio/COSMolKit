# RDKit Python EnumerateStereoisomers Full Source Port Plan

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
- Every real task step must be immediately preceded by reading `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
- The reading step must explicitly reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
- `Implement`, `Port`, `Modify`, `Update`, and `Fix` steps must produce a concrete artifact.
- `Audit` steps must produce a written gap report and must not replace implementation steps.
- If a step adds or updates tests, the next real task after the required reading step must run the most specific relevant test command for those tests.
- Do not defer tests added for one behavior to a final whole-plan validation step.
- Final whole-plan validation is still required when the plan changes code, but it does not replace immediate targeted validation after test-writing steps.
- If the plan violates this contract, regenerate the plan before doing any work.
- Copying C++ or Python comments, adding a dispatch stub, or adding placeholder branches is not a completed `Port` step.
- Do not use "smallest subpart", skeleton code, dispatch-only code, placeholder code, TODO-only code, or partial porting as completion.
- If a step is too large to complete fully, regenerate the plan by splitting that step into smaller full-port steps.
- Each split step must still represent a complete source-backed behavior, not a placeholder.
- A `Port` step is complete only when the selected source behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Exact Parity Boundary

- The public parity mainline is RDKit 2026.03.1 `rdkit.Chem.EnumerateStereoisomers.EnumerateStereoisomers()` and `GetStereoisomerCount()` from `third_party/rdkit/rdkit/Chem/EnumerateStereoisomers.py`, using the pinned reference identity in `testdata/reference/rdkit.json`.
- The potential-stereo engine is the C++ implementation bound as `Chem.FindPotentialStereo`: `third_party/rdkit/Code/GraphMol/Wrap/ChiralityOps.cpp`, `third_party/rdkit/Code/GraphMol/Chirality.h`, `third_party/rdkit/Code/GraphMol/FindStereo.cpp`, and the active cross-file helpers in `third_party/rdkit/Code/GraphMol/Chirality.cpp`.
- The Python path calls `FindPotentialStereo(mol, cleanIt=false, flagPossible=true)`; the exact internal engine must nevertheless preserve and test all four `cleanIt` and `flagPossible` combinations because they share one source implementation and cleanup is required to prove its mutation boundary.
- Atom, double-bond, and enhanced-stereo-group setter behavior follows `_AtomFlipper`, `_BondFlipper`, and `_StereoGroupFlipper` in the Python file, with the source-identical setter bodies in `third_party/rdkit/Code/GraphMol/EnumerateStereoisomers/Flippers.cpp` used as secondary line anchors.
- The 2025 C++ `StereoisomerEnumerator` is not the public parity target because its cleanup mode, conformer wedging, defaults, random iteration, RNG, uniqueness key, and atropisomer support differ from the Python API.
- `AtropisomerFlipper` is outside the Python enumeration boundary; already represented atropisomer state must remain valid and preserved where the Python source preserves it, but it must not be silently enumerated by the Python-parity API.
- RDKit's `tryEmbedding` filter is source-defined heuristic behavior and must be reproduced exactly as part of parity; no additional feasibility, small-ring, geometric-validity, or chemistry filter may be inserted before or after it.
- Official RDKit is a test-data oracle only and must not become a production dependency, FFI path, external-command fallback, or runtime subprocess.

## Architectural Convergence Requirements

- One typed potential-stereo engine owns atom, bond, ring, symmetry, cleanup, descriptor, controlling-atom, and specified-state decisions.
- One private flipper representation consumes typed potential-stereo records and enhanced stereo groups; no enumerator-specific candidate detector may survive.
- One lazy enumeration engine owns exhaustive configurations, random configurations, deduplication, embedding, output finalization, and count semantics.
- Project-native public APIs use `Molecule` value semantics and do not mutate the source molecule; internal cleanup and flippers operate only on an owned private enumeration workspace.
- Internal workspace mutation must not expose raw mutable `Molecule` storage or bypass the operation policy for any public mutation API.
- `StereoInfo.centeredOn` and `Atom::NOATOM` must become typed Rust state, not ambiguous `usize` values or magic integers.
- Existing exact shared implementations such as `count_swaps_to_interconvert`, `rank_fragment_atoms`, `symmetrize_sssr`, bridgehead detection, ring-stereo helpers, SMILES writing, hydrogen addition, and distance geometry must be reused after source-equivalence review instead of copied.
- Different-source behavior may retain a narrow adapter, but source-identical chemistry must have one implementation and one set of behavioral markers.
- The current `find_tetrahedral_centers`, `find_stereo_bonds`, `is_valid_double_bond_config`, eager `generate_combinations`, xorshift RNG, 20-center caps, and three duplicated enumeration loops are not valid completion paths.
- Compatibility wrappers may remain only when they delegate to the single core and have an intentional documented semantic boundary; undocumented parallel chemistry implementations must be removed.

## Required Observable State

- Potential-stereo oracle rows compare status, ordered record count, record order, stereo type, typed center, specified state, descriptor, permutation, ordered controlling atoms including missing-atom positions, atom and bond stereo state after analysis, atom and bond directions, stereo atoms, unknown-stereo state, computed-property membership, chiral ranks, ring-stereo state, and source-molecule preservation.
- Enumeration oracle rows compare options, selected configuration sequence, flipper kind and order, theoretical count, emitted count, emitted order, canonical isomeric SMILES, full atom and bond state, stereo groups, molecule properties including `_MolFileChiralFlag`, computed-property clearing, conformer count, embedding status, and coordinates for successful embedding.
- Discrete fields compare exactly; embedding coordinates use only the tolerance already justified for the pinned distance-geometry implementation.
- Every source failure or unsupported state maps to a structured project error; no panic, silent no-op, empty placeholder result, row filtering, expected mismatch, or heuristic fallback may support the parity claim.

## Required Test Domains

- Fixed source cases include every relevant `FindPotentialStereo` section in `third_party/rdkit/Code/GraphMol/catch_chirality.cpp`, every enumeration case in `third_party/rdkit/rdkit/Chem/UnitTestMol3D.py`, and the referenced RDKit fixture files with recorded provenance.
- Focused coverage includes tetrahedral atoms, implicit and explicit H, P, As, S, Se, ring N, bridgeheads, isotopes, non-tetrahedral centers, double bonds, cumulenes, para-stereo, fused and spiro rings, cages, unknown stereo, wiggly bonds, `STEREOANY`, enhanced stereo groups, no-center molecules, and already assigned centers.
- Option coverage includes every combination that changes control flow for `onlyUnassigned`, `onlyStereoGroups`, `unique`, `tryEmbedding`, `maxIsomers`, default deterministic sampling, integer seeds, and compatible custom Python random sources.
- Corpus coverage includes the repository small corpus, the 5,000-row corpus, and a bounded source-meaningful ChEMBL 37 phase that prevents combinatorial output explosion without filtering chemistry rows.
- Composition coverage includes parsed and constructed molecules, scalar and parallel calls, repeated calls, count followed by enumeration, SMILES roundtrip, existing conformers, source clone isolation, and cache-warm versus cache-cold execution.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit the pinned RDKit Python and C++ call graph, transitive helpers, current Rust implementations, public surfaces, operation boundaries, shared source ports, tests, fixtures, and exact observable-state gaps and write `dev/gap_reports/rdkit_python_enumerate_stereoisomers_source_audit.md`.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Correct premature behavioral markers and public comments in the current stereo enumeration path so no known heuristic, omitted branch, or non-corresponding source block remains marked behaviorally exact before replacement.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add the focused potential-stereo and stereoisomer-enumeration fixture corpus with source provenance, option matrices, branch rationale, referenced RDKit fixtures, and explicit in-scope and out-of-scope state.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Add the pinned-RDKit expected-data generator, schema, manifest inputs, and `generate_all.py` profile entries for exact potential-stereo and stereoisomer-enumeration oracle records.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Run the focused Python generator and schema tests for potential-stereo and stereoisomer-enumeration expected data.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Generate and validate every focused potential-stereo and stereoisomer-enumeration expected-data row with the pinned RDKit environment.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Implement the project-native typed stereo information model, typed atom-or-bond center, missing controlling-atom representation, specified state, descriptor, permutation, analysis options, result type, and structured error boundary.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Add focused model tests for every typed stereo information variant, equality and ordering behavior, missing controlling atoms, index validity, and serialization policy where applicable.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Run the exact typed stereo information model tests under strict operation contracts.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Port the active `getAtomNonzeroDegree`, protium-neighbor, tetrahedral-center, non-tetrahedral-center, and atom-potential-stereo source behavior by consolidating exact existing helpers instead of duplicating them.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Add focused atom-eligibility tests for all source element, degree, hydrogen, valence, charge, hybridization, conjugation, ring-N, bridgehead, query, and non-tetrahedral branches.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Run the exact atom-potential-stereo eligibility tests under strict operation contracts.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Port atom `getStereoInfo` with adjacency-order capture, sorted controlling atoms, swap parity, unknown stereo, tetrahedral descriptors, non-tetrahedral type inference, and chiral permutation semantics.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Add focused atom `StereoInfo` tests covering specified, unspecified, unknown, reordered neighbors, isotope state, missing ligands, non-tetrahedral permutations, and exact oracle record order.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Run the exact atom `StereoInfo` tests under strict operation contracts.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port bond-potential-stereo eligibility and bond `getStereoInfo` with degree and ring-size gates, ordered endpoint neighbors, missing-atom padding, unknown and wiggly state, E/Z translation, stereo-atom normalization, and represented atropisomer records.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Add focused bond `StereoInfo` tests covering invalid degrees, chain and ring double bonds, stereo-atom order, `STEREOANY`, `EITHERDOUBLE`, squiggle state, implicit-H padding, cumulene boundaries, and represented atropisomer bonds.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run the exact bond-potential-stereo and bond `StereoInfo` tests under strict operation contracts.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Port atom and bond comparison symbols plus `initAtomInfo` and `initBondInfo` with exact property-cache preparation, known and possible bitsets, source decorations, cleanup mutations, and represented-atrop preservation.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Add focused initialization tests for every specified-state branch, `cleanIt` and `flagPossible` combination, property-cache state, source symbol bytes, and cleanup side effect.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Run the exact potential-stereo initialization tests under strict operation contracts.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Port `flagRingStereo` and `areStereobondControllingAtomsDupes` with exact symmetrized-ring order, opposite-atom rules, divisor-two and divisor-three cases, fused-edge traversal, possible-ring counts, and atropisomer tie breaking.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Add focused ring-transmission tests for even, odd, fused, spiro, cage, para-stereo, opposite-center, opposite-bond, and over-aggressive historical regression cases.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run the exact ring-stereo transmission and controlling-atom duplicate tests under strict operation contracts.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Port `updateAtoms` and `updateBonds` with exact canonical-rank duplicate elimination, descriptor inversion, fixed-center state, symbol evolution, possible-center removal, ring-count rollback, controlling-atom ordering, and repeat-round signaling.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Add focused iterative-update tests for atom-only, bond-only, atom-to-bond, bond-to-atom, recursive para-stereo, ring rollback, descriptor inversion, and convergence-order cases.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Run the exact iterative atom and bond update tests under strict operation contracts.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Port `cleanMolStereo`, `atomIsCandidateForRingStereochem`, and `findChiralAtomSpecialCases` by consolidating source-identical existing ring helpers and reproducing exact atom tags, bond state, direction cleanup, stereo atoms, and computed ring properties.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Add focused cleanup and special-case tests for invalid marked atoms and bonds, slash-direction cleanup, non-tetrahedral permutation zeroing, ring-stereo candidate caching, BFS relationship signs, and unrelated-state preservation.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run the exact cleanup and ring special-case tests under strict operation contracts.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Port `runCleanup` and `findPotentialStereo` with exact symmetrized-ring preparation, non-strict property-cache update, custom-symbol canonical ranking, fixed-point passes, cleanup retry, chiral-rank state, potential-stereo cache state, result ordering, and structured failures.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Add full-field golden tests for every focused `findPotentialStereo` row and all `cleanIt` and `flagPossible` combinations without filtering or accepted mismatch.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Run the complete focused `findPotentialStereo` golden suite in release mode with strict operation contracts.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Implement the project-native value-style potential-stereo analysis surface on an owned private workspace so read-only analysis preserves the source molecule and cleanup state is returned explicitly rather than hidden mutation.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Add public potential-stereo API tests for source preservation, cleaned-result ownership, COW isolation, cache-warm and cache-cold equivalence, typed output, structured errors, and operation-policy compliance.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run the public potential-stereo API tests in release mode with strict operation contracts.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Port the Python `_AtomFlipper`, `_BondFlipper`, `_StereoGroupFlipper`, and `_getFlippers` behavior into one private typed flipper core with exact setter direction, stereo-atom initialization, skip conditions, option filtering, group order, and no atropisomer enumeration.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Add focused flipper tests for unassigned and assigned atoms, double bonds with and without stereo atoms, missing controlling atoms, absolute and relative stereo groups, `onlyUnassigned`, `onlyStereoGroups`, repeated flips, and represented atropisomer preservation.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Run the exact private flipper and flipper-selection tests under strict operation contracts.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Port `_RangeBitsGenerator`, `_UniqueRandomBitsGenerator`, arbitrary-center count semantics, duplicate-configuration rejection, and finite exhaustion without eagerly materializing `2^N` configurations.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Add focused configuration-source tests for zero centers, bit order, more than machine-word centers, exact theoretical counts, duplicate random draws, finite exhaustion, and allocation complexity.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Run the exact exhaustive and unique-random configuration-source tests.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Port the Python default molecule-invariant seed, integer-seed behavior, `random.Random` bit sequence, and compatible custom `getrandbits` source boundary without retaining the current xorshift implementation.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Add exact random-sequence tests for default seeding, atom-order invariance, integer seeds, the RDKit deterministic custom-random fixture, large center counts, repeated calls, and parallel isolation.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Run the exact Python-compatible random source and seed tests.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Port the Python enumerator preprocessing and configuration application path with exact source cloning, CIP clearing, unknown and either-double direction clearing, chiral-flag assignment, flipper order, and source mutation order.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Add focused preprocessing and configuration-application tests for no-center molecules, fully assigned molecules, partial assignment, `STEREOANY`, either-double fixtures, enhanced stereo groups, pre-existing conformers, properties, and source preservation.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Run the exact enumerator preprocessing and configuration-application tests under strict operation contracts.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Port stereoisomer finalization with exact stereo-group removal, double-bond neighbor directions, computed-property clearing while preserving rings, stereochemistry assignment parameters, canonical isomeric SMILES uniqueness, duplicate continuation, and output property state.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Add focused finalization tests for canonical duplicates, meso cases, cumulenes, enhanced groups, CIP reset, `_ChiralityPossible`, ring-cache preservation, output ordering, and all full-state expected fields.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Run the exact stereoisomer finalization and uniqueness tests in release mode with strict operation contracts.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Port `tryEmbedding` composition with exact hydrogen addition, configuration-derived 31-bit seed, one-conformer embedding call, failure filtering, heavy-atom coordinate copy, conformer append behavior, and finite configuration exhaustion.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Add focused embedding tests for successful and failed configurations, max-isomer counting after failures, no-infinite-loop exhaustion, existing conformers, coordinate copying, deterministic seeds, and source preservation.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Run the exact focused `tryEmbedding` stereoisomer tests in release mode with strict operation contracts.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Implement the single lazy Rust stereoisomer iterator, exact option defaults, exact count API, structured iterator errors, project-native `Molecule` methods, and deliberate facade re-exports over the unified potential-stereo and flipper core.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Add complete Rust integration tests for the pinned RDKit enumeration cases, option matrix, iterator laziness, no-center single result, arbitrary center counts, output order, repeated calls, composition calls, and scalar versus parallel behavior.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Run the complete focused Rust stereoisomer-enumeration integration suite in release mode with strict operation contracts.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Remove or reduce the legacy tetrahedral-only, double-bond-only, and combined enumeration paths to deliberate thin delegates while deleting heuristic candidate detection, heuristic ring filtering, eager combination allocation, center caps, xorshift sampling, duplicate loops, and non-source strategies.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Add architecture guard tests proving one potential-stereo owner, one flipper owner, one configuration engine, no direct public mutable-storage access, no retired helper use, and no duplicate source function bodies.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Run the stereo enumeration architecture, operation-policy, source-marker, and duplicate-source-port guards.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Implement the project-native Python molecule methods, options object, typed potential-stereo records, lazy iterator, count conversion, seeded random behavior, compatible custom random source, and structured exception mapping over the Rust core.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Build the Python extension in release mode with project-level `uv` and `maturin` tooling.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Add Python API tests for project-native naming, option defaults, iterator behavior, custom random sources, exact RDKit outputs, source preservation, errors, repr and typing, and mixed calls with SMILES, hydrogen, conformer, and descriptor APIs.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Run the focused Python potential-stereo and stereoisomer-enumeration tests against the release-built extension.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Regenerate `python/cosmolkit.pyi` from the implemented Python API without manual stub edits.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Run `basedpyright` over the focused Python tests and examples using the regenerated stubs.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Add repository small-corpus and 5,000-row potential-stereo and bounded stereoisomer-enumeration parity profiles with exact intermediate and final comparison fields.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Add parity-profile tests for manifest identity, complete row coverage, option coverage, comparison counts, restart safety, and rejection of filtered, partial, or mismatching output.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Run the small-corpus and 5,000-row parity-profile workflow tests.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Run every small-corpus and 5,000-row potential-stereo and bounded stereoisomer-enumeration parity profile and require zero unexplained mismatch.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Add a source-meaningful bounded potential-stereo and stereoisomer-enumeration phase to the repository-owned ChEMBL 37 parity workflow without excluding chemistry rows or permitting combinatorial output growth.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Add ChEMBL workflow tests for stereo phase identity, complete shard coverage, option and branch profiles, comparison totals, restart safety, and rejection of filtered or partial results.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Run the focused ChEMBL parity workflow tests.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Run the complete ChEMBL 37 potential-stereo and bounded stereoisomer-enumeration phase with every prepared shard and profile and require zero unexplained mismatch.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Run deterministic scalar, lazy-prefix, full-exhaustive, bounded-random, embedding, and parallel benchmarks against pinned RDKit while recording candidate discovery, configuration generation, finalization, and total costs.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Audit source-line closure, behavioral markers, algorithmic complexity, allocation shape, single-core ownership, operation-policy compliance, exact parity evidence, corpus completeness, and benchmark results and write `dev/gap_reports/rdkit_python_enumerate_stereoisomers_full_port_validation.md`.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Add a dedicated stereoisomer-enumeration `FeatureSpec` and update support metadata to `supported_with_rdkit_parity` only when the validation report records zero unexplained mismatch and every in-scope completion criterion is satisfied.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Update current Rust docs, Python docs, examples, README feature summaries, `VALIDATION.md`, `dev/porting_inventory.md`, and support tables with one consistent project-native API and final verified corpus evidence without editing historical plans.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Update the pending release section of `CHANGELOG.md` once with the completed potential-stereo and stereoisomer-enumeration port and final exact comparison totals.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Run `cargo fmt --all -- --check`.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Run `cargo clippy --workspace --all-targets --features cosmolkit-core/op-contracts-strict`.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Run the complete Python test suite against the release-built extension.
Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Run the Sphinx documentation build and doctest suite with warnings treated as errors.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Run repository source-marker, duplicate-source-port, generated-stub, support-matrix, expected-data-manifest, temporary-artifact, and documentation-link guards.
Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [x]: Move this fully completed plan to `dev/archive/plans/` without rewriting its checked execution history.
Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [x]: Update the active and archived plan indexes to reflect the completed plan's final location.
