# DG Bounds RDKit Source-Closure Plan

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
- Every real task step must be immediately preceded by reading `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md`.
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
- The baseline for this plan is the DG bounds generation scope defined by direct source comparison against `third_party/rdkit/Code/DistGeom/BoundsMatrix.h`, `third_party/rdkit/Code/DistGeom/TriangleSmooth.cpp`, `third_party/rdkit/Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp`, and the bounds-wrapper portion of `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp` through `getMolBoundsMatrix(...)`.
- The target definition of near-`100%` DG bounds source closure for this plan is: every bounds-generation function in that baseline has a source-backed Rust implementation or an explicit documented out-of-scope note, every remaining `RDKit✔️❌` or `RDKit❗✔️` marker in the DG bounds call chain is either upgraded with evidence or justified by a documented scope note, every new targeted test slice passes, and the final source audit finds no silent direct gap inside the selected bounds-generation scope.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit the current Rust DG bounds implementation against `BoundsMatrix.h`, `TriangleSmooth.cpp`, `BoundsMatrixBuilder.cpp`, and `rdDistGeom.cpp::getMolBoundsMatrix(...)` and write a fresh gap report to `dev/gap_reports/dg_bounds_source_closure_reaudit.md`.

Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Port `DistGeom::BoundsMatrix::setUpperBoundIfBetter(...)`, `DistGeom::BoundsMatrix::setLowerBoundIfBetter(...)`, and `DistGeom::BoundsMatrix::checkValid()` into `crates/cosmolkit-core/src/chemistry/distgeom.rs` with exact source anchors from `BoundsMatrix.h`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for `setUpperBoundIfBetter(...)`, `setLowerBoundIfBetter(...)`, and `checkValid()`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom bounds_matrix check_valid if_better`.

Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port the `DistGeom::triangleSmoothBounds(BoundsMatPtr, double)` forwarding overload into `crates/cosmolkit-core/src/chemistry/distgeom.rs` and align the pointer and pointer-like smoother entrypoints with `TriangleSmooth.cpp`.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for the new `triangleSmoothBounds(BoundsMatPtr, double)` forwarding behavior.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom triangle_smooth`.

Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Port `UFF::Tools::addAtomChargeFlags(...)` exactly in `crates/cosmolkit-core/src/chemistry/distgeom.rs` by replacing the remaining summary-style placeholder lines with direct source-backed logic and copied C++ comments.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for `addAtomChargeFlags(...)` charge and lanthanide label behavior.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom atom_charge_flags`.

Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Port `UFF::Tools::getAtomLabel(...)` exactly in `crates/cosmolkit-core/src/chemistry/distgeom.rs` by replacing the remaining summary-style placeholder lines with direct source-backed logic and copied C++ comments.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for `getAtomLabel(...)` hybridization, aromaticity, and charge-flag label composition.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom atom_label`.

Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Port the helper currently described as `RDGeom::compute13Dist via law of cosines` in `crates/cosmolkit-core/src/chemistry/distgeom.rs` with an exact source frame from the RDKit geometry helper that defines the 1-3 distance computation used by `_set13BoundsHelper(...)`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for the `compute13Dist` helper behavior used by `_set13BoundsHelper(...)`.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom compute13`.

Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port the shortest-topological-distance helper used to build the `distMatrix` input for `set14Bounds(...)` and `set15Bounds(...)` with direct source anchors from the RDKit helper call chain that supplies those distances.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for the shortest-topological-distance helper behavior used by `set14Bounds(...)` and `set15Bounds(...)`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom topological_distance`.

Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Upgrade the source reproduction frames and behavior closure for `set12Bounds(...)` in `crates/cosmolkit-core/src/chemistry/distgeom.rs` so that the remaining `RDKit❗✔️` marker lines are either upgraded with evidence or replaced by exact source-backed logic.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for any new `set12Bounds(...)` branches completed in Step 41.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set12_bounds`.

Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Upgrade the source reproduction frames and behavior closure for `set13Bounds(...)` in `crates/cosmolkit-core/src/chemistry/distgeom.rs` so that the remaining `RDKit✔️❌` marker lines are either upgraded with evidence or replaced by exact source-backed logic.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for any new `set13Bounds(...)` branches completed in Step 47.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set13_bounds`.

Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Upgrade the source reproduction frames and behavior closure for `_record14Path(...)`, `_setInRing14Bounds(...)`, `_setTwoInSameRing14Bounds(...)`, `_setTwoInDiffRing14Bounds(...)`, and `_setShareRingBond14Bounds(...)` in `crates/cosmolkit-core/src/chemistry/distgeom.rs`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for any new ring-path 1-4 branches completed in Step 53.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom in_ring14 same_ring14 diff_ring14 share_ring14`.

Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Upgrade the source reproduction frames and behavior closure for `_checkH2NX3H1OX2(...)`, `_checkNhChChNh(...)`, `_checkAmideEster14(...)`, `_checkMacrocycleAllInSameRingAmideEster14(...)`, `_isCarbonyl(...)`, and `_checkAmideEster15(...)` in `crates/cosmolkit-core/src/chemistry/distgeom.rs`.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for any new 1-4 and 1-5 classification-helper branches completed in Step 59.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom amide ester carbonyl macrocycle`.

Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Upgrade the source reproduction frames and behavior closure for `_setMacrocycleTwoInSameRing14Bounds(...)`, `_setMacrocycleAllInSameRing14Bounds(...)`, and `_setChain14Bounds(...)` in `crates/cosmolkit-core/src/chemistry/distgeom.rs`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for any new macrocycle and chain 1-4 branches completed in Step 65.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom macrocycle14 chain14 force_trans_amides`.

Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Upgrade the source reproduction frames and behavior closure for `set14Bounds(...)`, `_set15BoundsHelper(...)`, and `set15Bounds(...)` in `crates/cosmolkit-core/src/chemistry/distgeom.rs`.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for any new `set14Bounds(...)`, `_set15BoundsHelper(...)`, and `set15Bounds(...)` branches completed in Step 71.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set14_bounds set15_bounds set15_bounds_helper compute15`.

Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Upgrade the source reproduction frames and behavior closure for `collectBondsAndAngles(...)`, both `setTopolBounds(...)` overloads, and `RDKit::getMolBoundsMatrix(...)` in `crates/cosmolkit-core/src/chemistry/distgeom.rs`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Add or update targeted tests in `crates/cosmolkit-core/src/chemistry/distgeom/tests.rs` for any new wrapper and overload branches completed in Step 77.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom collect_bonds_and_angles set_topol_bounds dg_bounds_matrix`.

Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Re-audit `crates/cosmolkit-core/src/chemistry/distgeom.rs` and write `dev/gap_reports/dg_bounds_source_closure_final_audit.md` with a function-by-function verdict covering every remaining `RDKit✔️❌`, `RDKit❗✔️`, and `RDKit✔️❗` marker in the DG bounds call chain.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Update `dev/porting_inventory.md`, `dev/gap_reports/dg_bounds_remaining_source_scan.md`, and `crates/cosmolkit-core/src/support.rs` to reflect the exact final DG bounds source-closure state without overstating parity.

Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Run `cargo fmt --all`.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom -- --nocapture`.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --lib`.
