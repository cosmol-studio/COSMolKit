# RDKit 2D Coordinate Generation Full-Port Checklist

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
- The full-port baseline for this checklist is the RDKit 2D depiction source surface in:
  `third_party/rdkit/Code/GraphMol/Depictor/RDDepictor.h`,
  `third_party/rdkit/Code/GraphMol/Depictor/RDDepictor.cpp`,
  `third_party/rdkit/Code/GraphMol/Depictor/DepictUtils.h`,
  `third_party/rdkit/Code/GraphMol/Depictor/DepictUtils.cpp`,
  `third_party/rdkit/Code/GraphMol/Depictor/EmbeddedFrag.h`,
  `third_party/rdkit/Code/GraphMol/Depictor/EmbeddedFrag.cpp`,
  `third_party/rdkit/Code/GraphMol/Depictor/Templates.h`,
  `third_party/rdkit/Code/GraphMol/Depictor/Templates.cpp`,
  and `third_party/rdkit/Code/GraphMol/Depictor/Basement/Depictor.cpp`.
- The active Rust implementation surface for this checklist is expected to include, as needed,
  `crates/cosmolkit-core/src/chemistry/coordinates.rs`,
  `crates/cosmolkit-core/src/operations/ops.rs`,
  `crates/cosmolkit-core/src/model/molecule.rs`,
  `crates/cosmolkit-core/src/support.rs`,
  `crates/cosmolkit-core/src/properties/batch.rs`,
  `crates/cosmolkit-core/src/io/molblock.rs`,
  and targeted tests adjacent to the coordinate implementation.
- The target definition of `~100% RDKit 2D coordinate generation port coverage` for this checklist is:
  every baseline direct function and direct helper that contributes to 2D coordinate generation semantics is either fully ported with copied source anchors and targeted tests, or is explicitly documented in the final audit as a narrow non-Rust-runtime compatibility wrapper with no missing chemistry behavior.
- “Fully ported” in this checklist means:
  no simplified stub remains,
  no TODO-only completion remains,
  no missing helper in the active call chain remains,
  and the relevant targeted tests for the just-ported behavior pass before continuing.
- If a later audit step finds a new direct helper or direct function gap inside this baseline, the checklist must be regenerated before claiming near-100% coverage.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit the full RDKit 2D depiction baseline against the current Rust code and write `dev/gap_reports/coordinate_2d_remaining_source_scan.md` with a direct function-by-function inventory and remaining gap list.

Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Port `RDDepict::Compute2DCoordParameters`, `RDDepict::compute2DCoords(RDKit::ROMol &, const Compute2DCoordParameters &)`, and `RDDepict::compute2DCoords(RDKit::ROMol &, const RDGeom::INT_POINT2D_MAP *, bool, bool, unsigned int, unsigned int, int, bool, bool, bool)` into explicit Rust parameter and entrypoint APIs without placeholder branches.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add or update targeted tests for the Rust `Compute2DCoordParameters` mapping and both `compute2DCoords`-equivalent entrypoints, including default-argument behavior and non-default parameter forwarding.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates compute2d_params -- --nocapture`.

Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port `RDDepict::preferCoordGen`, `Compute2DCoordParameters::forceRDKit`, and `RDKit::Add2DCoordsToMol(ROMol &, bool useDLL)` semantics into Rust with source-backed behavior and explicit non-CoordGen runtime handling instead of implicit omission.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Add or update targeted tests for `preferCoordGen` and `forceRDKit` routing semantics and the `Add2DCoordsToMol`-equivalent wrapper behavior.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates coordgen_force_rdkit add2dcoords_to_mol -- --nocapture`.

Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.2 [x]: Port a dedicated Rust parser/helper for RDKit ring-template SMARTS/CX lines that splits the SMARTS body from the CX coordinate block and rejects line forms that are not supported by RDKit `SmartsToMol(line)` template loading semantics.
Step 17.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.4 [x]: Add or update targeted tests for template-line splitting and rejection behavior, including preserved CX coordinate text and source-backed failure for malformed or unsupported line forms.
Step 17.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.6 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates template_line_parse -- --nocapture`.
Step 17.7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.8 [x]: Port the minimal Rust runtime graph model needed to hold a ring-system template parsed from the SMARTS body, including atom queries, bond queries, explicit template-bond endpoints, and ring-closure-expanded connectivity.
Step 17.9 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.10 [x]: Add or update targeted tests for the template graph model, including ring-closure-expanded bond connectivity and bond-count preservation from RDKit-style SMARTS template rows.
Step 17.11 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.12 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates template_graph_model -- --nocapture`.
Step 17.13 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.14 [x]: Port the minimal Rust runtime coordinate/ring wrapper needed to combine the parsed template graph with preserved 2D coordinates and ring-info initialization sufficient for later `assertValidTemplate(...)` checks.
Step 17.15 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.16 [x]: Add or update targeted tests for the combined template runtime model, including SMARTS/CX coordinate retention and initialized ring-bond counts for valid ring-system template inputs.
Step 17.17 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17.18 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates template_runtime_model -- --nocapture`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port `RDDepict::CoordinateTemplates::assertValidTemplate(...)` exactly on top of the new template runtime model, including missing-coordinates, 3D-coordinate, multi-fragment, single-atom, and non-ring-bond rejection behavior.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Add or update targeted tests for `assertValidTemplate(...)`, including each explicit RDKit rejection message class.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates template_validation -- --nocapture`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25.2 [x]: Port the file-open and line-iteration shell of `RDDepict::CoordinateTemplates::loadTemplatesFromPath(...)` into a dedicated Rust loader helper with copied source anchors and the exact RDKit file-open failure message.
Step 25.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25.4 [x]: Port the per-line SMARTS parse failure path of `loadTemplatesFromPath(...)`, including the exact `"Could not load templates from {templatePath}: Invalid smarts"` error contract on the first invalid line.
Step 25.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25.6 [x]: Port the validated insertion path of `loadTemplatesFromPath(...)`, including runtime-model construction, ring-info initialization parity, `assertValidTemplate(...)` call chaining, and per-size append order within the destination bucket.
Step 25.7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25.8 [x]: Add or update targeted tests for template-path loading, including file-open failure, invalid-smarts failure, template-validation failure forwarding, and duplicate-template preservation order within the returned size bucket.
Step 25.9 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25.10 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates template_loading -- --nocapture`.

Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35.2 [x]: Port the Rust template-registry storage shell corresponding to `CoordinateTemplates::getRingSystemTemplates()`, `hasTemplateOfSize(...)`, and `getMatchingTemplates(...)`, including singleton-style process-global access and size-bucket queries.
Step 35.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35.4 [x]: Port `CoordinateTemplates::loadDefaultTemplates()` exactly on top of the Rust template runtime model and embedded RDKit `TEMPLATE_SMARTS`, including clear-and-reload behavior and invalid-default-entry skip semantics.
Step 35.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35.6 [x]: Port `CoordinateTemplates::setRingSystemTemplates(...)` and `RDDepict::setRingSystemTemplates(...)`, including load-into-temporary-then-replace semantics so current templates remain unchanged on load failure.
Step 35.7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35.8 [x]: Port `CoordinateTemplates::addRingSystemTemplates(...)`, `RDDepict::addRingSystemTemplates(...)`, and `RDDepict::loadDefaultRingSystemTemplates()`, including per-bucket front insertion order and default-template reset behavior.
Step 35.9 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35.10 [x]: Add or update targeted tests for set/add/default ring template APIs, including replacement semantics, failure-preserves-current-registry semantics, front-insertion order for added templates, and default-template reload behavior.
Step 35.11 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35.12 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates ring_system_templates -- --nocapture`.

Step 40.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.2 [x]: Port `RDDepict::embedRing(...)` and `RDDepict::transformPoints(...)` exactly into source-backed Rust helpers with copied source blocks.
Step 40.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.4 [x]: Add or update targeted tests for `embedRing` and `transformPoints`, including polygon radius and affine transform cases.
Step 40.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.6 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates embed_ring transform_points -- --nocapture`.
Step 40.7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.8 [x]: Port `RDDepict::computeBisectPoint(...)`, `RDDepict::reflectPoint(...)`, and `RDDepict::reflectPoints(...)` exactly into source-backed Rust helpers with copied source blocks.
Step 40.9 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.10 [x]: Add or update targeted tests for `computeBisectPoint`, `reflectPoint`, and `reflectPoints`, including obtuse-angle bisector and mirror cases.
Step 40.11 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.12 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates bisect reflect -- --nocapture`.

Step 46.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.2 [x]: Port `RDDepict::rankAtomsByRank(...)` exactly into the Rust helper surface, including RDKit CIP-rank-driven ordering, depict-rank fallback composition, stable ascending/descending ordering, and copied source anchors.
Step 46.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.4 [x]: Add or update targeted tests for `rankAtomsByRank(...)`, including CIP-rank precedence, depict-rank fallback ordering, and descending-order behavior.
Step 46.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.6 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates rank_atoms -- --nocapture`.
Step 46.7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.8 [x]: Port `RDDepict::setNbrOrder(...)` exactly into the Rust helper surface, including reference-neighbor injection, rank-based sorting, and the degree-4 opposite-side reorder behavior.
Step 46.9 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.10 [x]: Add or update targeted tests for `setNbrOrder(...)`, including the no-reference path and the reference-neighbor reorder path for degree-4 centers.
Step 46.11 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.12 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates set_nbr_order -- --nocapture`.
Step 46.13 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.14 [x]: Port `RDDepict::getNbrAtomAndBondIds(...)` and `RDDepict::findBondsPairsToPermuteDeg4(...)` exactly into the Rust helper surface with copied source anchors and degree-4 orthogonality checks.
Step 46.15 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.16 [x]: Add or update targeted tests for neighbor-id extraction and degree-4 bond-pair permutation selection, including perpendicular-neighbor and opposite-neighbor layouts.
Step 46.17 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 46.18 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates nbr_atom_bond_ids permute_deg4 -- --nocapture`.

Step 40.13 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.14 [x]: Port `RDDepict::pickFirstRingToEmbed(...)` exactly into the Rust fused-ring helper surface, including substituent counting and largest-ring tie breaking.
Step 40.15 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.16 [x]: Add or update targeted tests for `pickFirstRingToEmbed(...)`, including lower-substituent preference and larger-ring tie breaking.
Step 40.17 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.18 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates pick_first_ring -- --nocapture`.
Step 40.19 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.20 [x]: Port `RDDepict::findCoreRings(...)` exactly into the Rust fused-ring helper surface, including iterative one-atom / one-bond side-ring removal semantics.
Step 40.21 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.22 [x]: Add or update targeted tests for `findCoreRings(...)`, including removal of singly fused side rings and preservation of multi-ring cores.
Step 40.23 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.24 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates find_core_rings -- --nocapture`.
Step 40.25 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.26 [x]: Port `RDDepict::findNextRingToEmbed(...)` exactly into the Rust fused-ring helper surface, including the two-common-atoms early return path and common-atom chain rotation behavior.
Step 40.27 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.28 [x]: Add or update targeted tests for `findNextRingToEmbed(...)`, including two-common-atoms preference and reordered common-atom chain output for bridged cases.
Step 40.29 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 40.30 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates find_next_ring -- --nocapture`.

Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.2 [x]: Port `RDDepict::getAllRotatableBonds(...)` and `RDDepict::getRotatableBonds(...)` exactly into Rust helper functions, removing the current local closure-only duplication.
Step 47.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.4 [x]: Add or update targeted tests for `getAllRotatableBonds(...)` and `getRotatableBonds(...)`, including shortest-path trimming and non-ring/non-stereo filtering.
Step 47.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.6 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates rotatable_bonds -- --nocapture`.
Step 47.7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.8 [x]: Port `RDDepict::hasTerminalRGroupOrQueryHydrogen(...)` and `RDDepict::prepareTemplateForRGroups(...)` exactly into the Rust template/query helper surface, including former-index metadata propagation and terminal-removal semantics.
Step 47.9 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.10 [x]: Add or update targeted tests for terminal R-group/query-hydrogen detection and reduced-template preparation semantics.
Step 47.11 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.12 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates rgroups -- --nocapture`.
Step 47.13 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.14 [ ]: Port `RDDepict::reducedToFullMatches(...)` exactly into the Rust query-match helper surface, including in-place former-index restoration and unmatched-neighbor fill-in semantics.
Step 47.15 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.16 [ ]: Add or update targeted tests for `reducedToFullMatches(...)`, including restored former indices and appended neighbor matches.
Step 47.17 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.18 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates reduced_matches -- --nocapture`.
Step 47.19 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.20 [ ]: Port `RDDepict::invertWedgingIfMolHasFlipped(...)` exactly into Rust helper form and wire it to the active depiction transform surface where supported.
Step 47.21 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.22 [ ]: Add or update targeted tests for wedging inversion threshold behavior.
Step 47.23 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47.24 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates invert_wedging -- --nocapture`.

Step 52 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [ ]: Port `EmbeddedFrag::EmbeddedFrag(unsigned int aid, const ROMol *)`, `EmbeddedFrag::EmbeddedFrag(const ROMol &, const RDGeom::INT_POINT2D_MAP &)`, `EmbeddedFrag::EmbeddedFrag(const ROMol &, const INT_VECT &)`, `EmbeddedFrag::EmbeddedFrag(const Bond *)`, `EmbeddedFrag::computeNbrsAndAng(...)`, `EmbeddedFrag::findNumNeigh(...)`, `EmbeddedFrag::updateNewNeighs(...)`, `EmbeddedFrag::setupNewNeighs()`, `EmbeddedFrag::findNeighbor(...)`, and `EmbeddedFrag::setupAttachmentPoints()` exactly into Rust fragment state machinery.
Step 54 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [ ]: Add or update targeted tests for `EmbeddedFrag` constructors, neighbor/angle setup, attachment-point setup, and double-bond fragment seeding.
Step 56 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates embedded_frag_ctor compute_nbrs attachment_points -- --nocapture`.

Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.1 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.2 [x]: Port `checkStereoChemistry(...)` and `EmbeddedFrag::matchToTemplate(...)` fully, including exact-template-size filtering, ring-count/bond-count prechecks, degree-count screening, substructure-match extraction, stereo validation, and coordinate transfer into fragment state.
Step 59.3 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.4 [x]: Add or update targeted tests for stereo-preserving ring-template matching, including cis/trans acceptance, mismatch rejection, and template prefilter behavior.
Step 59.5 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.6 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates match_to_template stereo_template_match -- --nocapture`.
Step 59.7 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.8 [x]: Port `mirrorTransRingAtoms(...)`, `EmbeddedFrag::initFromRingCoords(...)`, and `EmbeddedFrag::mergeRing(...)` fully into the Rust fragment/ring helper surface, including trans-ring reflection and pinned-atom angle/nbr updates.
Step 59.9 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.10 [x]: Add or update targeted tests for trans-ring mirroring, ring-fragment initialization, and `mergeRing(...)` pinned-atom neighbor replacement behavior.
Step 59.11 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.12 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates mirror_trans_ring init_ring_coords merge_ring -- --nocapture`.
Step 59.13 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.14 [x]: Port `EmbeddedFrag::embedFusedRings(...)` fully on top of the Rust template registry and fragment helper surface, including whole-system template attempts, core-ring template fallback, first-ring seeding, iterative ring attachment, and done-ring bookkeeping.
Step 59.15 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.16 [x]: Add or update targeted tests for fused-ring embedding with and without templates, including whole-system-template hit, core-ring fallback, and first-ring seeding behavior.
Step 59.17 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59.18 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates embed_fused_rings -- --nocapture`.

Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Port `EmbeddedFrag::computeOneAtomTrans(...)`, `EmbeddedFrag::computeTwoAtomTrans(...)`, `EmbeddedFrag::Reflect(...)`, `EmbeddedFrag::reflectIfNecessaryCisTrans(...)`, `EmbeddedFrag::reflectIfNecessaryThirdPt(...)`, `EmbeddedFrag::reflectIfNecessaryDensity(...)`, `EmbeddedFrag::findCommonAtoms(...)`, `EmbeddedFrag::mergeNoCommon(...)`, `EmbeddedFrag::mergeWithCommon(...)`, `EmbeddedFrag::mergeFragsWithComm(...)`, `EmbeddedFrag::Transform(...)`, `EmbeddedFrag::computeBox()`, and `EmbeddedFrag::canonicalizeOrientation()` exactly, removing the current one-atom and third-point merge omissions.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Add or update targeted tests for one-atom and two-atom transforms, cis/trans reflection, third-point reflection, density-based reflection, fragment merging, box computation, and canonical orientation.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates one_atom_trans two_atom_trans reflect_if_necessary merge_frags canonicalize_orientation -- --nocapture`.

Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Port `EmbeddedFrag::addNonRingAtom(...)`, `EmbeddedFrag::addAtomToAtomWithAng(...)`, `EmbeddedFrag::addAtomToAtomWithNoAng(...)`, and `EmbeddedFrag::expandEfrag(...)` completely into the Rust acyclic expansion path.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Add or update targeted tests for non-ring atom addition, angle-based expansion, no-angle expansion, and `expandEfrag(...)` growth order.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates add_non_ring_atom add_atom_with_ang add_atom_with_no_ang expand_efrag -- --nocapture`.

Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Port `_recurseAtomOneSide(...)`, `_crossVal(...)`, `_pairDIICompAscending(...)`, `_findClosestPair(...)`, `EmbeddedFrag::computeDistMat(...)`, `EmbeddedFrag::mimicDistMatAndDensityCostFunc(...)`, `EmbeddedFrag::permuteBonds(...)`, `EmbeddedFrag::randomSampleFlipsAndPermutations(...)`, `EmbeddedFrag::findCollisions(...)`, `EmbeddedFrag::totalDensity()`, `_recurseDegTwoRingAtoms(...)`, `_anyNonRingBonds(...)`, `EmbeddedFrag::flipAboutBond(...)`, `_findDeg1Neighbor(...)`, `_findClosestNeighbor(...)`, `EmbeddedFrag::openAngles(...)`, `EmbeddedFrag::removeCollisionsBondFlip()`, `EmbeddedFrag::removeCollisionsOpenAngles()`, and `EmbeddedFrag::removeCollisionsShortenBonds()` completely, replacing the current adapted collision and sampling shortcuts with source-backed behavior.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Add or update targeted tests for mimic-distance cost, bond permutation, random sampling, collision detection, density scoring, bond flipping, open-angle repair, and bond-shortening collision cleanup.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates mimic_distmat permute_bonds random_sample find_collisions remove_collisions -- --nocapture`.

Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Port `DepictorLocal::getRankedAtomNeighbors(...)`, `DepictorLocal::embedSquarePlanar(...)`, `DepictorLocal::embedTBP(...)`, `DepictorLocal::embedOctahedral(...)`, and `DepictorLocal::embedNontetrahedralStereo(...)` completely and wire them into the actual 2D coordinate generation pipeline instead of leaving them disconnected.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Add or update targeted tests for non-tetrahedral stereo seeding and pipeline integration, including square-planar, trigonal-bipyramidal, and octahedral cases.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates nontetrahedral_stereo square_planar tbp octahedral -- --nocapture`.

Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Port `DepictorLocal::embedFusedSystems(...)`, `DepictorLocal::embedCisTransSystems(...)`, `DepictorLocal::getNonEmbeddedAtoms(...)`, `DepictorLocal::_findLargestFrag(...)`, `DepictorLocal::_shiftCoords(...)`, `copySign(...)`, and `ThetaBin` completely, removing the current fused-ring and cis/trans stubs.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Add or update targeted tests for fused-system seeding, cis/trans system embedding, non-embedded atom selection, largest-fragment selection, coordinate shifting, and `ThetaBin` binning behavior.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates embed_fused_systems embed_cis_trans get_non_embedded_atoms find_largest_frag shift_coords thetabin -- --nocapture`.

Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Port `computeInitialCoords(...)` and `copyCoordinate(...)` exactly, including ring discovery, stereochemistry setup, prespecified-coordinate fragment creation, fused-ring embedding, non-tetrahedral stereo embedding, cis/trans embedding, fragment expansion order, and conformer-copy semantics.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Add or update targeted tests for `computeInitialCoords(...)` and `copyCoordinate(...)`, including `coordMap`, `preSpec`, fused-ring-first, cis/trans-first, and conformer-clearing behavior.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates compute_initial_coords copy_coordinate -- --nocapture`.

Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Port the full top-level `compute2DCoords(...)` behavior, including random sampling vs collision-bond-flip selection, open-angle and shorten-bond cleanup, canonical orientation, `_shiftCoords(...)`, single-atom `coordMap` translation, and exact default-parameter control flow.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Add or update targeted tests for the full top-level `compute2DCoords(...)` path, including `coordMap`, `canonOrient`, `clearConfs`, `nFlipsPerSample`, `nSamples`, `sampleSeed`, `permuteDeg4Nodes`, `forceRDKit`, and `useRingTemplates`.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates compute2d_full coord_map canon_orient clear_confs ring_templates -- --nocapture`.

Step 106 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [ ]: Port `compute2DCoordsMimicDistMat(...)` exactly, including distance-matrix weighting, random sampling call chain, canonical orientation, and copied-conformer semantics.
Step 108 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [ ]: Add or update targeted tests for `compute2DCoordsMimicDistMat(...)`, including negative-entry ignore behavior, weight handling, and deterministic seed behavior.
Step 110 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates mimic_distmat_api -- --nocapture`.

Step 112 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [ ]: Port `removeAllConformersButOne(...)` and `generateDepictionMatching2DStructure(ROMol &, const ROMol &, const MatchVectType &, int, const ConstrainedDepictionParams &)` completely, including hard-constraint handling, `acceptFailure`, `alignOnly`, `allowRGroups`, and wedging-adjustment semantics.
Step 114 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [ ]: Add or update targeted tests for constrained 2D depiction with explicit `MatchVectType`, including failure acceptance, align-only mode, R-group handling, and molblock-wedging adjustment.
Step 116 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates constrained_2d_matchvect remove_all_conformers -- --nocapture`.

Step 118 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Port `generateDepictionMatching2DStructure(ROMol &, const ROMol &, const MatchVectType &, int, bool)`, `generateDepictionMatching2DStructure(ROMol &, const ROMol &, int, const ROMol *, const ConstrainedDepictionParams &)`, `generateDepictionMatching2DStructure(ROMol &, const ROMol &, int, const ROMol *, bool, bool, bool)`, and `generateDepictionMatching3DStructure(...)` completely, including overload-specific matching and reference-pattern semantics.
Step 120 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Add or update targeted tests for the remaining constrained-2D overloads and `generateDepictionMatching3DStructure(...)`, including reference-pattern and accept-failure branches.
Step 122 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates constrained_2d_overloads constrained_3d -- --nocapture`.

Step 124 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Port `straightenDepiction(...)` and `normalizeDepiction(...)` exactly, including angle histogram/binning, minimum-rotation handling, canonicalization control, and bond-length normalization scaling.
Step 126 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Add or update targeted tests for `straightenDepiction(...)` and `normalizeDepiction(...)`, including `minimizeRotation`, `canonicalize`, and explicit `scaleFactor` behavior.
Step 128 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict coordinates straighten_depiction normalize_depiction -- --nocapture`.

Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Port the public Rust exposure for the completed RDKit 2D depiction surface through `Molecule` APIs, registered operations, batch helpers, and existing `molblock`/drawing callers so they use the completed parameterized implementation instead of the current narrow default-only path.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Add or update targeted tests for public `with_2d_coordinates` operation behavior, value semantics, batch behavior, and the call sites that previously depended on the old default-only coordinate path.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict with_2d_coordinates batch.with_2d_coordinates molblock draw -- --nocapture`.

Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Audit the complete 2D depiction baseline against RDKit source again after the implementation steps and rewrite `dev/gap_reports/coordinate_2d_remaining_source_scan.md` with the exact remaining gap state and direct helper inventory.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Update `dev/porting_inventory.md` and `crates/cosmolkit-core/src/support.rs` to reflect the exact 2D coordinate generation port-completion state without unsupported parity overclaims.

Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Run `cargo fmt --all`.
Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.

Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Audit the final 2D coordinate generation baseline against RDKit source and write `dev/gap_reports/coordinate_2d_frozen_coverage_final_audit.md` stating whether the 2D depiction port is now at or near `100%` coverage for the selected baseline and listing any residual non-chemistry wrapper exclusions explicitly.
