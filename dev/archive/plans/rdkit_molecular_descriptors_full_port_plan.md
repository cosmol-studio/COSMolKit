# RDKit Molecular Descriptors Full Port Plan

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
- Do not use "smallest subpart", skeleton code, dispatch-only code, placeholder code, TODO-only code, or partial porting as completion.
- If a step is too large to complete fully, regenerate the plan by splitting that step into smaller full-port steps.
- Each split step must still represent a complete source-backed behavior, not a placeholder.
- A `Port` step is complete only when the selected RDKit behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Scope

- RDKit oracle version: the pinned `third_party/rdkit` tree used by the repository parity framework.
- Rust test surface: `crates/cosmolkit-core/tests/rdkit_molecular_descriptor_parity.rs`.
- Rust implementation surface: `crates/cosmolkit-core/src/properties/descriptors.rs`.
- Golden generator surface: `tests/scripts/gen_rdkit_molecular_descriptors_golden.py`.
- Golden artifacts: `tests/golden/smiles_small/molecular_descriptors.jsonl` and `tests/golden/smiles_5000/molecular_descriptors.jsonl`.
- Source-level target functions:
  - `third_party/rdkit/Code/GraphMol/Descriptors/MolDescriptors.cpp::calcAMW`
  - `third_party/rdkit/Code/GraphMol/Descriptors/MolDescriptors.cpp::calcExactMW`
  - `third_party/rdkit/Code/GraphMol/Descriptors/MolDescriptors.cpp::calcMolFormula`
  - `third_party/rdkit/Code/GraphMol/MolProps.cpp::getAvgMolWt`
  - `third_party/rdkit/Code/GraphMol/MolProps.cpp::getExactMolWt`
  - `third_party/rdkit/Code/GraphMol/MolProps.cpp::getMolFormula`
  - `third_party/rdkit/Code/GraphMol/MolProps.cpp::HillCompare`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcNumHBD`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcNumHBA`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcFractionCSP3`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcNumAromaticRings`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcNumRotatableBonds`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Crippen.cpp::getCrippenAtomContribs`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Crippen.cpp::calcCrippenDescriptors`
  - `third_party/rdkit/Code/GraphMol/Descriptors/Crippen.cpp::CrippenParamCollection::CrippenParamCollection`
  - `third_party/rdkit/Code/GraphMol/Descriptors/MolSurf.cpp::getTPSAAtomContribs`
  - `third_party/rdkit/Code/GraphMol/Descriptors/MolSurf.cpp::calcTPSA`
  - `third_party/rdkit/rdkit/Chem/QED.py::properties`
  - `third_party/rdkit/rdkit/Chem/QED.py::ads`
  - `third_party/rdkit/rdkit/Chem/QED.py::qed`

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit `crates/cosmolkit-core/src/properties/descriptors.rs`, `crates/cosmolkit-core/tests/rdkit_molecular_descriptor_parity.rs`, and `tests/scripts/gen_rdkit_molecular_descriptors_golden.py` against the Scope section and write `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Modify `tests/scripts/gen_all_rdkit_goldens.py` and `tests/README.md` if Step 3 finds any descriptor golden registration or condition-naming gap.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Run `.venv/bin/python tests/scripts/gen_all_rdkit_goldens.py --python .venv/bin/python --profile smiles_small --suite default --clean --jobs 4`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptor_golden_has_one_record_per_smiles -- --nocapture`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port `third_party/rdkit/Code/GraphMol/MolProps.cpp::getAvgMolWt` and `third_party/rdkit/Code/GraphMol/Descriptors/MolDescriptors.cpp::calcAMW` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_mol_wt`.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Port `third_party/rdkit/Code/GraphMol/MolProps.cpp::getExactMolWt` and `third_party/rdkit/Code/GraphMol/Descriptors/MolDescriptors.cpp::calcExactMW` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_exact_mol_wt`.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port `third_party/rdkit/Code/GraphMol/MolProps.cpp::HillCompare`, `third_party/rdkit/Code/GraphMol/MolProps.cpp::getMolFormula`, and `third_party/rdkit/Code/GraphMol/Descriptors/MolDescriptors.cpp::calcMolFormula` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_mol_formula`.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Audit current SMARTS parsing and substructure matching support for the `Lipinski.cpp::SMARTSCOUNTFUNC` HBD/HBA patterns and write the result into `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Port the complete source-backed behavior required by `Lipinski.cpp::SMARTSCOUNTFUNC(NumHBD, ...)` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_num_hbd` or into a private helper called only by that function.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Audit and, if necessary, source-upgrade the exact recursive-SMARTS parsing and substructure behavior required by `Lipinski.cpp::SMARTSCOUNTFUNC(NumHBA, ...)`, writing the result into `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 29A [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29B [x]: Source-upgrade the exact recursive-SMARTS root matching and bond-query behavior required by `Lipinski.cpp::SMARTSCOUNTFUNC(NumHBA, ...)`, with copied RDKit source anchors and focused tests.
Step 29C [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29D [x]: Port the complete source-backed behavior required by `Lipinski.cpp::SMARTSCOUNTFUNC(NumHBA, ...)` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_num_hba` or into a private helper called only by that function.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcFractionCSP3` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_fraction_csp3`.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcNumAromaticRings` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_num_aromatic_rings`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcNumRotatableBonds` `NonStrict` and `Strict` branches into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_num_rotatable_bonds`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/Lipinski.cpp::calcNumRotatableBonds` `StrictLinkages` branch into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_num_rotatable_bonds`.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Audit Crippen SMARTS parameter support by comparing `third_party/rdkit/Code/GraphMol/Descriptors/Crippen.cpp::defaultParamData`, `CrippenParamCollection::CrippenParamCollection`, and current COSMolKit SMARTS/substructure behavior and write the result into `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/Crippen.cpp::CrippenParamCollection::CrippenParamCollection` into a private static parameter table helper for `crates/cosmolkit-core/src/properties/descriptors.rs::calc_crippen_descriptors`.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/Crippen.cpp::getCrippenAtomContribs` into a private helper called by `crates/cosmolkit-core/src/properties/descriptors.rs::calc_crippen_descriptors`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/Crippen.cpp::calcCrippenDescriptors` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_crippen_descriptors`.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/MolSurf.cpp::getTPSAAtomContribs` into a private helper called by `crates/cosmolkit-core/src/properties/descriptors.rs::calc_tpsa`.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Port `third_party/rdkit/Code/GraphMol/Descriptors/MolSurf.cpp::calcTPSA` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_tpsa`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Audit QED dependencies by comparing `third_party/rdkit/rdkit/Chem/QED.py::properties`, `QED.py::ads`, `QED.py::qed`, and the descriptor functions already ported in this plan, then write the result into `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Port `third_party/rdkit/rdkit/Chem/QED.py::ads` and the QED constant tables into private helpers used by `crates/cosmolkit-core/src/properties/descriptors.rs::calc_qed`.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Port `third_party/rdkit/rdkit/Chem/QED.py::properties` into a private helper called by `crates/cosmolkit-core/src/properties/descriptors.rs::calc_qed`.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Port `third_party/rdkit/rdkit/Chem/QED.py::qed` into `crates/cosmolkit-core/src/properties/descriptors.rs::calc_qed`.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity molecular_descriptors_match_rdkit_golden_for_supported_properties -- --nocapture`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run `COSMOLKIT_PARITY_PROFILE=smiles_5000 cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity -- --nocapture`.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Run `cargo fmt --all`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_molecular_descriptor_parity -- --nocapture`.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Audit the post-Step-81 source markers in `crates/cosmolkit-core/src/properties/descriptors.rs`, `crates/cosmolkit-core/src/search/substruct.rs`, `tests/scripts/gen_all_rdkit_goldens.py`, and `tests/README.md`, and update `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md` with every remaining descriptor/QED dependency gap, including exact file/function/line anchors for `RDKit❌❌`, `RDKit❗*`, and `RDKit✔️❌` markers.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Modify `tests/scripts/gen_all_rdkit_goldens.py`, `tests/README.md`, and the descriptor/QED golden generator surface so focused `Chem.DeleteSubstructs` parity is a named, reusable golden surface with generation conditions encoded in the filename or profile path, without regenerating unrelated golden files for this step.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Generate the focused `Chem.DeleteSubstructs` golden files for both `smiles_small` and `smiles_5000` through the unified `delete-substructs` generator suite, and verify the generated `delete_substructs_onlyfrags_chirality.jsonl` records cover `onlyFrags=false`, `onlyFrags=true`, `useChirality=false`, and `useChirality=true`.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Port `third_party/rdkit/Code/GraphMol/ChemTransforms/ChemTransforms.cpp::deleteSubstructs` `onlyFrags=true` branch into `crates/cosmolkit-core/src/properties/descriptors.rs::rdkit_qed_delete_substructs`, including the copied RDKit source lines for `MolOps::getMolFrags`, sorted fragment/match comparison, `Union(mxi, delList, tmp)`, and whole-fragment deletion semantics.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Add or update the smallest stable Rust test that exercises `rdkit_qed_delete_substructs` against the focused `Chem.DeleteSubstructs` golden for `onlyFrags=false` and `onlyFrags=true`, comparing canonical isomeric SMILES, atom count, and bond count exactly, with no tolerance, pass-rate reporting, row skipping, or heuristic fallback.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict properties::descriptors::qed_tests::delete_substructs -- --nocapture`.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Audit `third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp`, `third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.h`, and `third_party/rdkit/Code/GraphMol/Substruct/vf2.hpp` for every `SubstructMatchParameters::useChirality` and `specifiedStereoQueryMatchesUnspecified` dependency required by `Chem.DeleteSubstructs(..., useChirality=true)`, and record the exact function-level port targets in `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Extend `crates/cosmolkit-core/src/search/substruct.rs::SubstructMatchParams` with source-backed RDKit fields required for this plan, starting with `use_chirality` and `specified_stereo_query_matches_unspecified`, and update all existing construction sites so default behavior remains RDKit-compatible.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Port the `useChirality` atom-label and bond-label functor branches from `SubstructMatch.cpp` into `crates/cosmolkit-core/src/search/substruct.rs`, including tetrahedral query/molecule chiral-tag presence checks and double-bond stereo presence checks with `specifiedStereoQueryMatchesUnspecified` semantics.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Port the `MolMatchFinalCheckFunctor` chirality final-check path from `SubstructMatch.cpp` into `crates/cosmolkit-core/src/search/substruct.rs`, including tetrahedral chirality consistency, double-bond stereo consistency, uniquification interaction, and structured unsupported errors for any RDKit stereochemistry state not modeled by COSMolKit.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Extend the focused `Chem.DeleteSubstructs` golden generator to include RDKit's chiral cases from `third_party/rdkit/Code/GraphMol/ChemTransforms/testChemTransforms.cpp`, covering `useChirality=false` and `useChirality=true` for matching and opposite-chirality SMARTS.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Run the focused `Chem.DeleteSubstructs` parity test after the `useChirality=true` golden extension and require exact canonical isomeric SMILES, atom count, and bond count parity for every row.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Audit the remaining Crippen cache/property markers in `crates/cosmolkit-core/src/properties/descriptors.rs::calc_crippen_descriptors`, `rdkit_parse_default_crippen_params`, and `rdkit_crippen_atom_contribs`, and decide per marker whether COSMolKit models the RDKit mutable property-cache state; write the decision and required artifact for each marker into `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: If Step 107 determines COSMolKit models Crippen descriptor property-cache state, port RDKit's `_crippenLogP`, `_crippenMR`, `_crippenLogPContribs`, and `_crippenMRContribs` read/write branches into the appropriate molecule/property-cache layer and add exact tests for cached-value reuse and cache population; otherwise replace the remaining `RDKit❌❌` markers with source-protocol-compliant modeled-input-state comments that explicitly exclude mutable RDKit property-cache input.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Audit the remaining TPSA cache/property markers in `crates/cosmolkit-core/src/properties/descriptors.rs::calc_tpsa` and `rdkit_tpsa_atom_contribs`, and decide per marker whether COSMolKit models RDKit mutable property-cache state and atom-contribution cache state; write the decision and required artifact for each marker into `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md`.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: If Step 111 determines COSMolKit models TPSA property-cache state, port RDKit's `_tpsa-*` and `_tpsaAtomContribs-*` read/write branches into the appropriate molecule/property-cache layer and add exact tests for `force=false`, `force=true`, and `includeSandP` cache-key behavior; otherwise replace the remaining `RDKit❌❌` markers with source-protocol-compliant modeled-input-state comments that explicitly exclude mutable RDKit property-cache input.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Audit `third_party/rdkit/rdkit/Chem/QED.py::properties` `mol is None` behavior against the COSMolKit public `calc_qed(&Molecule)` API and either add a public/API-level test proving the input state is unrepresentable or introduce a wrapper-level nullable input path that reproduces RDKit's `ValueError`, then update the source marker accordingly.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Update `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md` so it no longer contains stale statements such as "every descriptor function currently fails closed" after implementation, and so every remaining unsupported or excluded source behavior has an explicit blocker, modeled-state boundary, or follow-up step number.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity -- --nocapture`.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Run `COSMOLKIT_PARITY_PROFILE=smiles_5000 cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity -- --nocapture`.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Run `cargo fmt --all`.
Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_molecular_descriptor_parity -- --nocapture`.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Update `dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md` so every remaining descriptor/QED-adjacent marker after Step 127 has an explicit disposition with a blocker, modeled-state boundary, or concrete follow-up plan target, replacing vague "later" or "future" wording.
Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Run `rg -n "Required follow-up|later|future|not currently blocking|whole-surface|outside this descriptor plan" dev/gap_reports/rdkit_molecular_descriptors_source_inventory.md` and update the report if any remaining match is an unresolved completion claim.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity -- --nocapture`.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run `COSMOLKIT_PARITY_PROFILE=smiles_5000 cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_molecular_descriptor_parity -- --nocapture`.
Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
