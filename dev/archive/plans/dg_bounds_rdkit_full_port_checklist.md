# DG Bounds RDKit Full-Port Checklist

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
- The full-port baseline for this checklist is the function list in `dev/gap_reports/dg_bounds_remaining_source_scan.md`.
- The target definition of `100% DG bounds port coverage` for this checklist is: every baseline function is fully ported with direct RDKit source comparison, every corresponding local targeted test step has passed, and the final RDKit source rescan finds no remaining direct helper or direct function gap inside the DG bounds scope.
- "Fully ported" means the Rust implementation reproduces the current modeled RDKit behavior for the selected function and its direct helper lines, not a subset, not a simplified branch, not a reduced-mode implementation, and not a placeholder unsupported return.
- If a function depends on an unported direct helper inside the baseline, the step is not complete until that helper is also covered by an earlier completed step.
- If any final audit step finds a new direct RDKit helper or direct function gap inside this DG bounds scope, the checklist must be regenerated before claiming `100% port coverage`.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Read `dev/gap_reports/dg_bounds_remaining_source_scan.md`.

Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Port `initBoundsMat(DistGeom::BoundsMatrix *mmat, ...)` and `initBoundsMat(DistGeom::BoundsMatPtr mmat, ...)` into Rust with source-backed matrix initialization semantics and explicit callable Rust entrypoints.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add or update targeted tests for the Rust `initBoundsMat` equivalents, including default-min/default-max behavior and diagonal/non-diagonal matrix cases.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom init_bounds`.

Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port `_checkAndSetBounds` exactly, including `setIfBetter` semantics, and align the Rust bounds-matrix core with the RDKit lower/upper update rules.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Add or update targeted tests for `_checkAndSetBounds` behavior, including conservative updates, better-bound replacement, and no-op cases.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom check_and_set_bounds`.

Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Port `set12Bounds` completely, including the RDKit UFF-derived bond-rest-length call chain semantics required for source-equivalent 1-2 bounds.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Add or update targeted tests for `set12Bounds`, including conjugated 5-ring squish behavior and UFF-backed bond-length branches.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set12_bounds`.

Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Port `isHBondAcceptor`, `isHBondDonor`, `isHinHBondDonor`, and `setLowerBoundVDW` completely with RDKit-equivalent lower-bound scaling behavior.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Add or update targeted tests for the H-bond helpers and `setLowerBoundVDW`, including 1-5 scaling, 1-6 scaling, and full-VDW branches.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom hbond set_lower_bound_vdw`.

Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Port `isLargerSP2Atom`, `_setRingAngle`, and `_set13BoundsHelper` completely, including larger-SP2 tolerance widening and ring-angle mutation rules.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Add or update targeted tests for `isLargerSP2Atom`, `_setRingAngle`, and `_set13BoundsHelper`, including larger-SP2 and ring-size special cases.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set13_bounds_helper ring_angle sp2`.

Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port `set13Bounds` completely, including ring-first ordering, fused-ring handling, and all visited-bound semantics required by RDKit.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Add or update targeted tests for `set13Bounds`, including fused-ring and non-ring 1-3 path coverage.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set13_bounds`.

Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Port `_checkH2NX3H1OX2`, `_checkNhChChNh`, `_checkAmideEster14`, `_checkMacrocycleAllInSameRingAmideEster14`, `_isCarbonyl`, and `_checkAmideEster15` completely as standalone RDKit helper equivalents for the 1-4 and 1-5 call chains.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add or update targeted tests for the 1-4 and 1-5 classification helpers completed in Step 41.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom amide ester carbonyl macrocycle`.

Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Port `_record14Path`, `_setInRing14Bounds`, `_setTwoInSameRing14Bounds`, `_setTwoInDiffRing14Bounds`, and `_setShareRingBond14Bounds` completely with RDKit-equivalent path classification and bound-setting behavior.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Add or update targeted tests for the ring-based 1-4 helper behaviors completed in Step 47.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom in_ring14 same_ring14 diff_ring14 share_ring14`.

Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Port `_checkMacrocycleTwoInSameRingAmideEster14`, `_setMacrocycleTwoInSameRing14Bounds`, and `_setMacrocycleAllInSameRing14Bounds` completely with RDKit-equivalent macrocycle heuristics.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Add or update targeted tests for the macrocycle 1-4 helper behaviors completed in Step 53.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom macrocycle14`.

Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Port `_setChain14Bounds` completely, including double-bond stereo, sulfur special case, amide/ester forcing, and `forceTransAmides` behavior.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Add or update targeted tests for `_setChain14Bounds`, including chain double-bond, sulfur, amide, and force-trans-amide branches.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom chain14 force_trans_amides`.

Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Port `set14Bounds` completely, including ring-path ordering, macrocycle switch behavior, ring-bond-pair bookkeeping, done-path bookkeeping, and all helper dispatch behavior.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Add or update targeted tests for `set14Bounds`, including same-ring, two-ring, shared-middle-ring, chain, and macrocycle entrypoint coverage.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set14_bounds`.

Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Port `_compute15DistsCisCis`, `_compute15DistsCisTrans`, `_compute15DistsTransTrans`, and `_compute15DistsTransCis` completely with exact RDKit numerical behavior and argument ordering.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Add or update targeted tests for the four `_compute15...` helpers completed in Step 71.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom compute15`.

Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Port `_set15BoundsHelper` completely against the full RDKit path14 and cis/trans bookkeeping behavior.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Add or update targeted tests for `_set15BoundsHelper`, including visited-bound skips, unknown-path fallback, and set15Atoms bookkeeping.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set15_bounds_helper`.

Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Port `set15Bounds` completely against the RDKit entry loop and reverse-path dispatch behavior.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Add or update targeted tests for `set15Bounds`, including reverse-path handling and interaction with the completed 1-4 path set.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set15_bounds`.

Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Port `collectBondsAndAngles` completely, including the exact bond-pair and angle-record generation used by the RDKit overload.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Add or update targeted tests for `collectBondsAndAngles`, including triple-bond and consecutive-double-bond angle flags.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom collect_bonds_and_angles`.

Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Port the first `setTopolBounds` overload completely, including `set15bounds`, `scaleVDW`, `useMacrocycle14config`, `forceTransAmides`, `set14bounds`, and `set13bounds` control flow.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Add or update targeted tests for the first `setTopolBounds` overload, including disabled-stage and enabled-stage permutations.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set_topol_bounds`.

Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Port the second `setTopolBounds` overload completely, including `bonds` and `angles` output collection through the completed `collectBondsAndAngles` path.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Add or update targeted tests for the second `setTopolBounds` overload, including exact `bonds` and `angles` outputs.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom set_topol_bounds bonds angles`.

Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Replace the current `dg_bounds_matrix` wrapper with a source-backed wrapper that is explicitly built on the completed `initBoundsMat` and `setTopolBounds` path without retaining simplified behavior.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Add or update targeted tests for `Molecule::dg_bounds_matrix()` and the public wrapper behavior after the completed `setTopolBounds` port.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict dg_bounds_matrix`.

Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Port `DistGeom::BoundsMatrix` storage semantics from `third_party/rdkit/Code/DistGeom/BoundsMatrix.h` into the Rust DG bounds core, including exact upper-triangle/lower-triangle read-write behavior instead of the current split symmetric `lower`/`upper` representation.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Add or update targeted tests for the Rust `BoundsMatrix` representation, including exact lower-triangle export, upper-triangle export, `getLowerBound`, `getUpperBound`, `setLowerBound`, and `setUpperBound` behavior against RDKit semantics.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Run the most specific strict test command covering the new `BoundsMatrix` storage/export tests.

Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Re-port `_checkAndSetBounds(...)` on top of the new raw matrix layout if needed, so the Rust implementation still matches the RDKit `BoundsMatrix.h` access pattern after the storage rewrite.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Add or update targeted tests for `_checkAndSetBounds(...)` after the matrix-layout rewrite, including conservative updates and set-if-better overlap behavior through the new backing representation.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom check_and_set_bounds`.

Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Port `DistGeom::triangleSmoothBounds(...)` exactly from `third_party/rdkit/Code/DistGeom/TriangleSmooth.cpp`, including the RDKit lower-bound difference formulas, tolerance reconciliation branch, and failure condition.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Add or update targeted tests for the Rust triangle smoother, including exact lower-bound tightening, tolerance-based upper-bound adjustment, and inconsistent-bounds failure behavior.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Run the most specific strict test command covering the new triangle-smoothing tests.

Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Port the public DG bounds wrapper behavior from `third_party/rdkit/Code/GraphMol/DistGeomHelpers/Wrap/rdDistGeom.cpp::getMolBoundsMatrix(...)`, including default arguments, optional smoothing behavior, and raw matrix export shape.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Add or update targeted tests for `Molecule::dg_bounds_matrix()` and the Rust wrapper layer, including RDKit default `scaleVDW=false`, `set15bounds=true`, optional smoothing control, and exported lower-triangle semantics.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict dg_bounds_matrix`.

Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Update the legacy local `distgeom` unit tests that still assume pre-port behavior, including the current `set13Bounds requires ring info` failures, so they assert the source-backed DG bounds pipeline instead of the old simplified assumptions.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Run the most specific strict test command covering the repaired legacy `distgeom` unit tests that currently fail under the new source-backed path.

Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Audit the complete DG bounds call graph against RDKit source again after the matrix-layer, smoother-layer, and wrapper-layer work, and rewrite `dev/gap_reports/dg_bounds_remaining_source_scan.md` with a zero-gap result if and only if every baseline function is fully covered.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Update `dev/porting_inventory.md` and `crates/cosmolkit-core/src/support.rs` to reflect the exact DG bounds port-completion state without introducing parity claims.

Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Run `cargo fmt --all`.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.

Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Audit the DG bounds baseline against RDKit source one final time and write `dev/gap_reports/dg_bounds_frozen_coverage_final_audit.md` stating whether DG bounds port coverage is exactly `100%`.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Port `DGeomHelpers::_getAtomStereo(const Bond *bnd, unsigned int aid1, unsigned int aid4)` as an explicit direct RDKit baseline function with copied source lines in `crates/cosmolkit-core/src/chemistry/distgeom.rs`, and update any affected direct-helper marker framing so the DG bounds baseline function inventory is source-complete.
Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Add or update targeted tests for `_getAtomStereo(...)`, including stereo-preserving and stereo-flipping branches driven by `stereoAtoms`.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict distgeom get_atom_stereo -- --nocapture`.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Rewrite `dev/gap_reports/dg_bounds_remaining_source_scan.md` so the direct RDKit baseline function list explicitly includes `_getAtomStereo(...)`, removes the current direct-helper inventory omission, and refreshes the targeted `distgeom` test-count evidence to the current executed result.
Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Rewrite `dev/gap_reports/dg_bounds_source_closure_final_audit.md` and `dev/gap_reports/dg_bounds_frozen_coverage_final_audit.md` so the final `100%` conclusion is claimed only if the direct RDKit baseline inventory, source anchors, and current targeted validation evidence are fully consistent after Step 159.
