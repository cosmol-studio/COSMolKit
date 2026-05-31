# RDKit Conformer Generation Full Port Plan

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
- Every real task step must be immediately preceded by reading:
  `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
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

## Scope Contract

- The implementation target is full RDKit conformer generation coverage, not a subset and not a heuristic approximation.
- The source baseline is `third_party/rdkit/Code/DistGeom/*`, `third_party/rdkit/Code/GraphMol/DistGeomHelpers/*`, `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp`, and every force-field or CrystalFF function reached by RDKit `EmbedMolecule` and `EmbedMultipleConfs`.
- Already ported COSMolKit DG bounds and force-field functions may be reused only after a source-level reuse audit proves that the behavior needed by conformer generation is already covered.
- Any unsupported branch discovered during execution must be explicitly ported or explicitly represented as a structured unsupported error only when RDKit behavior itself cannot be represented by current COSMolKit state.
- The final public surface must include Rust and Python value-style APIs for single-conformer and multi-conformer generation, parameter presets, failure tracking, optional force-field cleanup, and example workflows without relying on RDKit interop.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit the complete RDKit conformer-generation source baseline and write `dev/gap_reports/rdkit_conformer_generation_source_inventory.md` listing every class, struct, enum, free function, method, data table, fixture, wrapper, and reachable force-field dependency.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit existing COSMolKit DG bounds and force-field code against the inventory and write `dev/gap_reports/rdkit_conformer_generation_reuse_map.md` with exact reused Rust symbols and remaining source gaps.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Port `DistGeom::ChiralSet` and related aliases from `third_party/rdkit/Code/DistGeom/ChiralSet.h` into the Rust conformer-generation module with copied source anchors and exhaustive unit tests.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom_chiral_set -- --nocapture`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port `DistGeom::calcChiralVolume` overloads from `ChiralViolationContribs.cpp` and `ChiralViolationContribs.h` with copied source anchors and exhaustive unit tests.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chiral_volume -- --nocapture`.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Port `DistGeom::ChiralViolationContribs::ChiralViolationContribs`, `addContrib`, `getEnergy`, and `getGrad` from `ChiralViolationContribs.cpp` with copied source anchors and exhaustive unit tests.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chiral_violation_contribs -- --nocapture`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port `DistGeom::DistViolationContribs::DistViolationContribs`, local `distance2`, local `distance`, `getEnergy`, and `getGrad` from `DistViolationContribs.cpp` with copied source anchors and exhaustive unit tests.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict dist_violation_contribs -- --nocapture`.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Port every reachable class and method from `third_party/rdkit/Code/DistGeom/FourthDimContribs.h` with copied source anchors and exhaustive unit tests.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict fourth_dim_contribs -- --nocapture`.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Port `DistGeom::pickRandomDistMat` seed overload and RNG overload from `DistGeomUtils.cpp` and `DistGeomUtils.h` with copied source anchors and exhaustive deterministic unit tests.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict pick_random_dist_mat -- --nocapture`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Port `DistGeom::computeInitialCoords` seed overload and RNG overload from `DistGeomUtils.cpp` and `DistGeomUtils.h` with copied source anchors and exhaustive eigenvalue, negative-eigenvalue, zero-eigenvalue, and seeded-random unit tests.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict compute_initial_coords -- --nocapture`.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port `DistGeom::computeRandomCoords` seed overload and RNG overload from `DistGeomUtils.cpp` and `DistGeomUtils.h` with copied source anchors and exhaustive range, seed, and empty-position unit tests.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict compute_random_coords -- --nocapture`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Port `DistGeom::constructForceField` from `DistGeomUtils.cpp` and `DistGeomUtils.h` with copied source anchors and exhaustive distance, chiral, fourth-dimension, fixed-point, extra-weight, and basin-threshold unit tests.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict construct_distgeom_forcefield -- --nocapture`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Port `DistGeom::addImproperTorsionTerms` from `DistGeomUtils.cpp` with copied source anchors and exhaustive unit tests.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom_improper_torsion_terms -- --nocapture`.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Port `DistGeom::addExperimentalTorsionTerms` from `DistGeomUtils.cpp` with copied source anchors and exhaustive unit tests.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom_experimental_torsion_terms -- --nocapture`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Port `DistGeom::add12Terms` from `DistGeomUtils.cpp` with copied source anchors and exhaustive unit tests.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom_add12_terms -- --nocapture`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Port `DistGeom::add13Terms` from `DistGeomUtils.cpp` with copied source anchors and exhaustive unit tests.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom_add13_terms -- --nocapture`.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Port `DistGeom::addLongRangeDistanceConstraints` from `DistGeomUtils.cpp` with copied source anchors and exhaustive unit tests.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom_long_range_distance_constraints -- --nocapture`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Port `DistGeom::construct3DForceField` overloads from `DistGeomUtils.cpp` and `DistGeomUtils.h` with copied source anchors and exhaustive unit tests.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict construct_3d_forcefield -- --nocapture`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Port `DistGeom::constructPlain3DForceField` from `DistGeomUtils.cpp` and `DistGeomUtils.h` with copied source anchors and exhaustive unit tests.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict construct_plain_3d_forcefield -- --nocapture`.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Port `DistGeom::construct3DImproperForceField` overloads from `DistGeomUtils.cpp` and `DistGeomUtils.h` with copied source anchors and exhaustive unit tests.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict construct_3d_improper_forcefield -- --nocapture`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Audit CrystalFF conformer-generation dependencies against `TorsionPreferences.cpp`, `TorsionAngleContribs.cpp`, `TorsionAngleM6.cpp`, and existing Rust CrystalFF code and write `dev/gap_reports/rdkit_conformer_generation_crystalff_reuse_audit.md`.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Port any missing `CrystalFF::calcTorsionEnergyM6` overloads from `TorsionAngleContribs.cpp` and `TorsionAngleM6.cpp` with copied source anchors and exhaustive unit tests.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict crystalff_torsion_energy_m6 -- --nocapture`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Port any missing `CrystalFF::TorsionAngleContribs::addContrib`, `getEnergy`, and `getGrad` behavior from `TorsionAngleContribs.cpp` with copied source anchors and exhaustive unit tests.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict crystalff_torsion_angle_contribs -- --nocapture`.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Port any missing `CrystalFF::TorsionAngleContribM6::getEnergy` and `getGrad` behavior from `TorsionAngleM6.cpp` with copied source anchors and exhaustive unit tests.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict crystalff_torsion_angle_m6 -- --nocapture`.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Port any missing `CrystalFF::getExperimentalTorsions` overloads from `TorsionPreferences.cpp` with copied source anchors and exhaustive torsion-library, small-ring, bridged-ring, macrocycle, ETversion, and SMARTS-match unit tests.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict crystalff_experimental_torsions -- --nocapture`.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Port `DGeomHelpers::EmbedFailureCauses` from `Embedder.h` with public Rust status mapping and exhaustive unit tests.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_failure_causes -- --nocapture`.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Port `DGeomHelpers::EmbedParameters` and its constructor defaults from `Embedder.h` with copied source anchors and exhaustive field-default unit tests.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_parameters_defaults -- --nocapture`.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Port global parameter presets `KDG`, `ETDG`, `ETDGv2`, `ETKDG`, `ETKDGv2`, `ETKDGv3`, and `srETKDGv3` from `Embedder.cpp` with copied source anchors and exhaustive preset field tests.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_parameter_presets -- --nocapture`.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Port `DGeomHelpers::updateEmbedParametersFromJSON` from `EmbedderUtils.cpp` with copied source anchors and exhaustive JSON update tests.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict update_embed_parameters_from_json -- --nocapture`.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Port `DGeomHelpers::embedParametersToJSON` from `EmbedderUtils.cpp` with copied source anchors and exhaustive JSON export tests.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_parameters_to_json -- --nocapture`.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Port local helper `DGeomHelpers::haveOppositeSign` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_have_opposite_sign -- --nocapture`.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Port failure-count mutex behavior `failmutex_get`, `failmutex_create`, and `GetFailMutex` from `Embedder.cpp` into Rust synchronization or a documented no-shared-mutex equivalent with copied source anchors and exhaustive concurrency tests.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_failure_mutex -- --nocapture`.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Port `_volumeTest` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_volume_test -- --nocapture`.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Port `_sameSide` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_same_side -- --nocapture`.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Port `_centerInVolume` atom-index overload and chiral-set overload from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_center_in_volume -- --nocapture`.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Port `_boundsFulfilled` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_bounds_fulfilled -- --nocapture`.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Port `generateInitialCoords` from `Embedder.cpp` with copied source anchors and exhaustive random-coordinate, distance-matrix, seed, negative-eigenvalue, and zero-eigenvalue unit tests.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_generate_initial_coords -- --nocapture`.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Port `firstMinimization` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_first_minimization -- --nocapture`.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Port `checkTetrahedralCenters` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_check_tetrahedral_centers -- --nocapture`.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Port `checkChiralCenters` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_check_chiral_centers -- --nocapture`.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Port `minimizeFourthDimension` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_minimize_fourth_dimension -- --nocapture`.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Port `minimizeWithExpTorsions` from `Embedder.cpp` with copied source anchors and exhaustive ETKDG, KDG, CPCI, bounds-scaling, and convergence unit tests.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_minimize_with_exp_torsions -- --nocapture`.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Port `doubleBondGeometryChecks` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_double_bond_geometry_checks -- --nocapture`.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Port `doubleBondStereoChecks` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_double_bond_stereo_checks -- --nocapture`.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Port `finalChiralChecks` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_final_chiral_checks -- --nocapture`.
Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Port `embedPoints` from `Embedder.cpp` with copied source anchors and exhaustive failure-cause, timeout, callback, retry-loop, force-field-cleanup, and ETKDG unit tests.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_embed_points -- --nocapture`.
Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [x]: Port `DGeomHelpers::findDoubleBonds` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_find_double_bonds -- --nocapture`.
Step 180 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 181 [x]: Port `DGeomHelpers::findChiralSets` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 182 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 183 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_find_chiral_sets -- --nocapture`.
Step 184 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 185 [x]: Port `DGeomHelpers::adjustBoundsMatFromCoordMap` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 186 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 187 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_adjust_bounds_mat_from_coord_map -- --nocapture`.
Step 188 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 189 [x]: Port `DGeomHelpers::initETKDG` from `Embedder.cpp` with copied source anchors and exhaustive basic-knowledge, torsion, improper, ring, macrocycle, and force-trans-amide unit tests.
Step 190 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 191 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_init_etkdg -- --nocapture`.
Step 192 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 193 [x]: Port `DGeomHelpers::setupInitialBoundsMatrix` from `Embedder.cpp` with copied source anchors and exhaustive smoothing, fallback, custom-bounds, macrocycle14, VDW-scaling, coord-map, and smoothing-failure unit tests.
Step 194 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 195 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_setup_initial_bounds_matrix -- --nocapture`.
Step 196 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 197 [x]: Port `_fillAtomPositions` from `Embedder.cpp` with copied source anchors and exhaustive unit tests.
Step 198 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 199 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_fill_atom_positions -- --nocapture`.
Step 200 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 201 [x]: Port `_isConfFarFromRest` from `Embedder.cpp` with copied source anchors and exhaustive RMS pruning, heavy-atom, symmetry, and terminal-group symmetrization unit tests.
Step 202 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 203 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_is_conf_far_from_rest -- --nocapture`.
Step 204 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 205 [x]: Port `multiplication_overflows_` from `Embedder.cpp` with copied source anchors and exhaustive overflow boundary unit tests.
Step 206 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 207 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_multiplication_overflows -- --nocapture`.
Step 208 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 209 [x]: Port `embedHelper_` from `Embedder.cpp` with copied source anchors and exhaustive multi-thread scheduling, sequential-random-seed, pruning, failure-tracking, callback, timeout, and conformer-ID unit tests.
Step 210 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 211 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_helper -- --nocapture`.
Step 212 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 213 [x]: Port `getMolSelfMatches` from `Embedder.cpp` with copied source anchors and exhaustive symmetry-pruning unit tests.
Step 214 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 215 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embedder_mol_self_matches -- --nocapture`.
Step 216 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 217 [x]: Port `DGeomHelpers::EmbedMultipleConfs(ROMol&, INT_VECT&, unsigned int, EmbedParameters&)` from `Embedder.cpp` and `Embedder.h` into value-style COSMolKit molecule APIs with copied source anchors and exhaustive unit tests.
Step 218 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 219 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_multiple_confs -- --nocapture`.
Step 220 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 221 [x]: Port inline `DGeomHelpers::EmbedMultipleConfs` return-vector overload from `Embedder.h` with copied source anchors and exhaustive unit tests.
Step 222 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 223 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_multiple_confs_return_vector -- --nocapture`.
Step 224 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 225 [x]: Port inline `DGeomHelpers::EmbedMolecule(ROMol&, EmbedParameters&)` from `Embedder.h` with copied source anchors and exhaustive unit tests.
Step 226 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 227 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_molecule_params -- --nocapture`.
Step 228 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 229 [x]: Port legacy inline `DGeomHelpers::EmbedMolecule` overloads from `Embedder.h` with copied source anchors and exhaustive unit tests.
Step 230 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 231 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_molecule_legacy_overloads -- --nocapture`.
Step 232 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 233 [x]: Port legacy inline `DGeomHelpers::EmbedMultipleConfs` overloads from `Embedder.h` with copied source anchors and exhaustive unit tests.
Step 234 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 235 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict embed_multiple_confs_legacy_overloads -- --nocapture`.
Step 236 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 237 [x]: Port wrapper behavior `rdDistGeom.cpp::EmbedMolecule` and `rdDistGeom.cpp::EmbedMolecule2` into the Python-facing Rust API design with copied source anchors and exhaustive wrapper unit tests.
Step 238 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 239 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rd_distgeom_embed_molecule_wrapper -- --nocapture`.
Step 240 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 241 [x]: Port wrapper behavior `rdDistGeom.cpp::EmbedMultipleConfs` and `rdDistGeom.cpp::EmbedMultipleConfs2` into the Python-facing Rust API design with copied source anchors and exhaustive wrapper unit tests.
Step 242 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 243 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rd_distgeom_embed_multiple_confs_wrapper -- --nocapture`.
Step 244 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 245 [x]: Port wrapper parameter factories `getKDG`, `getETDG`, `getETDGv2`, `getETKDG`, `getETKDGv2`, `getETKDGv3`, and `getsrETKDGv3` from `rdDistGeom.cpp` with copied source anchors and exhaustive unit tests.
Step 246 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 247 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rd_distgeom_parameter_factories -- --nocapture`.
Step 248 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 249 [x]: Port wrapper helper `getExpTorsHelper` from `rdDistGeom.cpp` with copied source anchors and exhaustive unit tests.
Step 250 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 251 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rd_distgeom_exp_tors_helper -- --nocapture`.
Step 252 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 253 [x]: Port wrapper helper `getExpTorsHelperWithParams` from `rdDistGeom.cpp` with copied source anchors and exhaustive unit tests.
Step 254 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 255 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rd_distgeom_exp_tors_helper_with_params -- --nocapture`.
Step 256 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 257 [x]: Port wrapper helper `embedParametersToJSONHelper` from `rdDistGeom.cpp` with copied source anchors and exhaustive unit tests.
Step 258 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 259 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rd_distgeom_embed_parameters_json_helper -- --nocapture`.
Step 260 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 261 [x]: Import RDKit conformer-generation fixtures from `third_party/rdkit/Code/DistGeom/testDistGeom.cpp`, `third_party/rdkit/Code/GraphMol/DistGeomHelpers/testDgeomHelpers.cpp`, `catch_tests.cpp`, `Wrap/testDistGeom.py`, and `test_data/` into COSMolKit parity fixtures with provenance metadata.
Step 262 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 263 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation_fixture_inventory -- --nocapture`.
Step 264 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 265 [x]: Add source-backed parity tests for single-conformer DG, KDG, ETDG, ETDGv2, ETKDG, ETKDGv2, ETKDGv3, and srETKDGv3 generation across the imported fixture corpus.
Step 266 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 267 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation_single_parity -- --nocapture`.
Step 268 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 269 [x]: Add source-backed parity tests for multi-conformer generation, pruning, symmetry pruning, sequential random seeds, failure tracking, timeout, and fragment embedding across the imported fixture corpus.
Step 270 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 271 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation_multi_parity -- --nocapture`.
Step 272 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 273 [x]: Add source-backed parity tests for coord-map, custom-bounds-matrix, CPCI, random-coordinate, and clear-conformer behavior across the imported fixture corpus. `ignoreSmoothingFailures` remains covered at the source-anchored branch level, but no synthetic coord-map parity fixture was added because RDKit coord maps are Euclidean distance maps and do not source-back an arbitrary smoothing-failure construction.
Step 274 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 275 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation_parameter_parity -- --nocapture`.
Step 276 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 277 [x]: Add source-backed parity tests for chirality, double-bond stereo, linear double bonds, final chiral bounds, final center-in-volume, macrocycle torsions, small-ring torsions, and trans-amide enforcement.
Step 278 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 279 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation_stereo_parity -- --nocapture`.
Step 280 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 281 [x]: Register conformer-generation support features, operation metadata, invariant requirements, and parity policy in the COSMolKit operation and support registries.
Step 282 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 283 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation_operation_registry -- --nocapture`.
Step 284 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 285 [x]: Expose Rust facade APIs for `EmbedParameters`, conformer-generation presets, `embed_molecule`, `embed_multiple_confs`, `with_3d_conformer`, and `with_3d_conformers`.
Step 286 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 287 [x]: Run `cargo test -p cosmolkit --features cosmolkit-core/op-contracts-strict conformer_generation_public_api -- --nocapture`.
Step 288 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 289 [x]: Expose Python APIs for conformer-generation parameters, presets, single-conformer generation, multi-conformer generation, failure tracking, and optional UFF/MMFF post-optimization.
Step 290 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 291 [x]: Run `cargo check -p cosmolkit-py`.
Step 292 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 293 [x]: Regenerate Python stubs with `cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen`.
Step 294 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 295 [x]: Run `.venv/bin/maturin develop --manifest-path python/Cargo.toml`.
Step 296 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 297 [x]: Add Python tests for single-conformer generation, multi-conformer generation, ETKDGv3 defaults, parameter presets, failure tracking, value semantics, and optional force-field post-optimization.
Step 298 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 299 [x]: Run `.venv/bin/pytest python/tests/test_python_api_sanity.py`.
Step 300 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 301 [x]: Add Rust and Python examples for native conformer generation, multi-conformer pruning, ETKDGv3 generation, and UFF/MMFF post-optimization without RDKit interop.
Step 302 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 303 [x]: Run `cargo run -p cosmolkit --example conformer_generation`.
Step 304 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 305 [x]: Run `.venv/bin/python python/examples/conformer_generation.py`.
Step 306 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 307 [x]: Update Rust README, root README, Python quickstart, Python molecule docs, and Python API reference for native conformer generation without overstating support beyond the completed source-backed surface.
Step 308 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 309 [x]: Run `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html`.
Step 310 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 311 [x]: Run `.venv/bin/basedpyright python/tests python/examples` (reported 0 errors and existing warning-only diagnostics; command exits nonzero under the current warning policy).
Step 312 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 313 [x]: Audit all conformer-generation Rust modules for remaining `RDKit❌❌`, `RDKit❗❗`, `RDKit❗✔️`, `RDKit✔️❌`, and `RDKit✔️❗` markers and write `dev/gap_reports/rdkit_conformer_generation_final_marker_audit.md`.
Step 314 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 315 [x]: Fix every actionable conformer-generation marker gap found in `dev/gap_reports/rdkit_conformer_generation_final_marker_audit.md` with source-backed code, tests, or an explicit structured unsupported artifact.
Step 316 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 317 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict conformer_generation -- --nocapture`.
Step 318 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 319 [x]: Run `cargo test --workspace --features cosmolkit-core/op-contracts-strict conformer_generation -- --nocapture`.
Step 320 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 321 [x]: Run `cargo test --workspace --features cosmolkit-core/op-contracts-strict`.
Step 322 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 323 [x]: Run `cargo fmt --all`.
Step 324 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 325 [x]: Update `dev/porting_inventory.md`, `crates/cosmolkit-core/src/support.rs`, and `dev/gap_reports/rdkit_conformer_generation_full_port_validation.md` with the exact completed conformer-generation support state and validation commands.
