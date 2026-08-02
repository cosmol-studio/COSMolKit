# SMILES RDKit Frozen Full-Port Checklist

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
- The frozen coverage baseline for this checklist is the function list in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`.
- The target definition of `100% coverage` for this checklist is: every frozen baseline function is fully ported with direct RDKit source comparison, every corresponding test step has passed, and the final RDKit source rescan finds no remaining direct helper or direct function gap inside the frozen scope.
- "Fully ported" means the Rust implementation reproduces the current modeled RDKit behavior for the selected function and its direct helper lines, not a subset, not a partial branch, not a reduced-mode implementation, and not a placeholder unsupported return.
- If a function depends on an unported direct helper inside the frozen scope, the step is not complete until that helper is also covered by an earlier completed step.
- If any final audit step finds a new direct RDKit helper or direct function gap inside the frozen scope, the checklist must be regenerated before claiming `100% coverage`.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Read `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`.

Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Port `ReportParseError`, `CleanupAfterParseError`, `ClearAtomChemicalProps`, and `AddFragToMol` with direct RDKit source reproduction and shared parser-path integration.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add or update targeted tests for `ReportParseError`, `CleanupAfterParseError`, `ClearAtomChemicalProps`, and `AddFragToMol`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict smiles_parse_ops`.

Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port `setDoubleBondNeighborDirections`, `setBondStereoFromDirections`, `clearDirFlags`, `clearAllBondDirFlags`, and `assignStereochemistryFrom3D` with shared core integration for SMILES and SDF call sites.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Add or update targeted tests for `setDoubleBondNeighborDirections`, `setBondStereoFromDirections`, `clearDirFlags`, `clearAllBondDirFlags`, and `assignStereochemistryFrom3D`.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict stereochemistry_from_3d`.

Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Port `SmilesLexer::new`, `SmilesLexer::next_token`, `parse_simple_atomd`, `parse_bracket_atomd`, `parse_charge_element`, `parse_h_element`, `parse_chiral_element`, `parse_element`, `parse_number`, and `parse_ring_number`.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Add or update targeted tests for `SmilesLexer::new`, `SmilesLexer::next_token`, `parse_simple_atomd`, `parse_bracket_atomd`, `parse_charge_element`, `parse_h_element`, `parse_chiral_element`, `parse_element`, `parse_number`, and `parse_ring_number`.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict parse_simple_atomd`.

Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Port `parse_mol`, `add_first_atom`, `add_atom_connected_to_active`, `add_explicit_bond_to_atom`, `add_single_bond_to_atom`, and `add_disconnected_atom`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Add or update targeted tests for `parse_mol`, `add_first_atom`, `add_atom_connected_to_active`, `add_explicit_bond_to_atom`, `add_single_bond_to_atom`, and `add_disconnected_atom`.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict add_first_atom`.

Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Port `add_branch_atom_connected_to_active`, `add_branch_explicit_bond`, `add_branch_single_bond`, `close_branch`, `add_ring_marker`, and `add_explicit_bond_ring_marker`.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Add or update targeted tests for `add_branch_atom_connected_to_active`, `add_branch_explicit_bond`, `add_branch_single_bond`, `close_branch`, `add_ring_marker`, and `add_explicit_bond_ring_marker`.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict close_branch`.

Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port `add_single_bond_ring_marker`, `close_ring_opening`, `check_ring_closure_branch_status`, `finish_parse`, and `close_mol_rings`.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Add or update targeted tests for `add_single_bond_ring_marker`, `close_ring_opening`, `check_ring_closure_branch_status`, `finish_parse`, and `close_mol_rings`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ring_closure`.

Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Port `check_chirality_specifications`, `set_unspecified_bond_types`, `adjust_atom_chirality_flags`, `cleanup_after_parsing`, `get_unspecified_bond_type_for_atoms`, and `nontetrahedral_chiral_permutation_for_probe`.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add or update targeted tests for `check_chirality_specifications`, `set_unspecified_bond_types`, `adjust_atom_chirality_flags`, `cleanup_after_parsing`, `get_unspecified_bond_type_for_atoms`, and `nontetrahedral_chiral_permutation_for_probe`.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict cleanup_after_parsing`.

Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Port `preprocess_smiles`, `handle_cx_part_and_name`, `to_mol`, and `mol_from_smiles`.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Add or update targeted tests for `preprocess_smiles`, `handle_cx_part_and_name`, `to_mol`, and `mol_from_smiles`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict handle_cx_part_and_name`.

Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Port `assign_double_bond_stereo_from_directions`, `assign_chiral_types_from_bond_dirs`, and `assign_chiral_types_from_3d`.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Add or update targeted tests for `assign_double_bond_stereo_from_directions`, `assign_chiral_types_from_bond_dirs`, and `assign_chiral_types_from_3d`.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict assign_chiral_types_from_3d`.

Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Port `parse_cx_extensions`, `parse_cx_coords`, `parse_cx_atom_labels`, `parse_cx_atom_values`, `parse_cx_atom_props`, `parse_cx_coordinate_bonds`, and `parse_cx_zero_bonds`.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Add or update targeted tests for `parse_cx_extensions`, `parse_cx_coords`, `parse_cx_atom_labels`, `parse_cx_atom_values`, `parse_cx_atom_props`, `parse_cx_coordinate_bonds`, and `parse_cx_zero_bonds`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict parse_cx_coords`.

Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Port `parse_cx_enhanced_stereo`, `parse_cx_unsaturation`, `parse_cx_ring_bonds`, `parse_cx_linknodes`, `parse_cx_data_sgroup_attr`, and `parse_cx_data_sgroup`.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Add or update targeted tests for `parse_cx_enhanced_stereo`, `parse_cx_unsaturation`, `parse_cx_ring_bonds`, `parse_cx_linknodes`, `parse_cx_data_sgroup_attr`, and `parse_cx_data_sgroup`.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict parse_cx_ring_bonds`.

Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Port `parse_cx_sgroup_hierarchy`, `parse_cx_polymer_sgroup`, and `parse_cx_substitution`.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Add or update targeted tests for `parse_cx_sgroup_hierarchy`, `parse_cx_polymer_sgroup`, and `parse_cx_substitution`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sgroup_hierarchy`.

Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Port `parse_cx_wedged_bonds`, `parse_cx_variable_attachments`, `set_cx_stereo_for_bond`, `parse_cx_doublebond_stereo`, `process_cx_radical_section`, `parse_cx_radicals`, and `remove_hs_update_explicit_count`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Add or update targeted tests for `parse_cx_wedged_bonds`, `parse_cx_variable_attachments`, `set_cx_stereo_for_bond`, `parse_cx_doublebond_stereo`, `process_cx_radical_section`, `parse_cx_radicals`, and `remove_hs_update_explicit_count`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict remove_hs_update_explicit_count`.

Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Port `mol_to_random_smiles_vect`, `mol_to_smiles_with_mode`, `prepare_plain_smiles_molecule`, and `collect_fragment_write_plans`.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Add or update targeted tests for `mol_to_random_smiles_vect`, `mol_to_smiles_with_mode`, `prepare_plain_smiles_molecule`, and `collect_fragment_write_plans`.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol_to_random_smiles_vect`.

Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Port `choose_fragment_start_atom`, `canonicalize_fragment_stack`, `write_mol_stack`, `mol_fragment_to_smiles`, and `assemble_fragment_smiles`.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Add or update targeted tests for `choose_fragment_start_atom`, `canonicalize_fragment_stack`, `write_mol_stack`, `mol_fragment_to_smiles`, and `assemble_fragment_smiles`.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol_fragment_to_smiles`.

Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Port `get_bond_smiles`, `get_molecule_bond_smiles`, `atom_needs_bracket`, and `write_cx_zero_bonds`.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Add or update targeted tests for `get_bond_smiles`, `get_molecule_bond_smiles`, `atom_needs_bracket`, and `write_cx_zero_bonds`.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict get_bond_smiles`.

Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Port `writer_parse_molfile_bond_cfg`, `write_cx_ringbond_cistrans_block`, `write_cx_linknodes_block`, `writer_parse_sgroup_crossings`, and `cleanup_stereo_groups_for_cx_smiles`.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Add or update targeted tests for `writer_parse_molfile_bond_cfg`, `write_cx_ringbond_cistrans_block`, `write_cx_linknodes_block`, `writer_parse_sgroup_crossings`, and `cleanup_stereo_groups_for_cx_smiles`.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict write_cx_ringbond_cistrans_block`.

Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Port `canonical_dfs_traversal`, `compute_writer_chiral_adjustments`, `insert_implicit_nontetrahedral_neighbors_for_writer`, and `chiral_atom_needs_tag_inversion_for_writer`.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Add or update targeted tests for `canonical_dfs_traversal`, `compute_writer_chiral_adjustments`, `insert_implicit_nontetrahedral_neighbors_for_writer`, and `chiral_atom_needs_tag_inversion_for_writer`.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict canonical_dfs_traversal`.

Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Port `canonicalize_double_bond_for_writer`, `handle_dir_conflicts_across_double_bond_for_writer`, `remove_unwanted_bond_dir_specs_for_writer`, `remove_redundant_bond_dir_specs_for_writer`, and `same_side_dirs_are_compatible_for_writer`.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Add or update targeted tests for `canonicalize_double_bond_for_writer`, `handle_dir_conflicts_across_double_bond_for_writer`, `remove_unwanted_bond_dir_specs_for_writer`, and `same_side_dirs_are_compatible_for_writer`.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict canonicalize_double_bond`.

Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Port `rank_fragment_atoms`, `init_fragment_canon_atoms`, `init_canon_atoms_stereo_group_assignment`, and `make_canon_bond_holder`.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Add or update targeted tests for `rank_fragment_atoms`, `init_fragment_canon_atoms`, `init_canon_atoms_stereo_group_assignment`, and `make_canon_bond_holder`.
Step 122 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rank_fragment_atoms`.

Step 124 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [x]: Port `rank_with_functor`, `refine_partitions`, `break_ties`, `hanoisort`, and `hanoi`.
Step 126 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [x]: Add or update targeted tests for `rank_with_functor`, `refine_partitions`, `break_ties`, `hanoisort`, and `hanoi`.
Step 128 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict hanoi`.

Step 130 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [x]: Port `special_chirality_atom_compare`, `special_symmetry_atom_compare`, `atom_compare_base`, and `compare_ring_atoms_concerning_num_neighbors`.
Step 132 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [x]: Add or update targeted tests for `special_chirality_atom_compare`, `special_symmetry_atom_compare`, `atom_compare_base`, and `compare_ring_atoms_concerning_num_neighbors`.
Step 134 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict compare_ring_atoms_concerning_num_neighbors`.

Step 136 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Port `symmetrize_sssr_with_options`, `find_ring_families`, and `find_sssr_internal`.
Step 138 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Add or update targeted tests for `symmetrize_sssr_with_options`, `find_ring_families`, and `find_sssr_internal`.
Step 140 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict symmetrize_sssr`.

Step 142 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [x]: Port `kekulize_assignment`, `kekulize_fragment_assignment`, and `kekulize_fused_assignment`.
Step 144 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [x]: Add or update targeted tests for `kekulize_assignment`, `kekulize_fragment_assignment`, and `kekulize_fused_assignment`.
Step 146 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_assignment`.

Step 148 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [x]: Port `mark_double_bond_candidates`, `kekulize_worker_matching`, `permute_dummies_and_match`, and `kekulize_if_possible_assignment`.
Step 150 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [x]: Add or update targeted tests for `mark_double_bond_candidates`, `kekulize_worker_matching`, `permute_dummies_and_match`, and `kekulize_if_possible_assignment`.
Step 152 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 153 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_if_possible`.

Step 154 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 155 [x]: Port `ValenceContext::new`, `ValenceContext::atom`, `ValenceContext::incident_bonds`, `periodic_table_row`, `required_valence_list`, and `rdkit_element_symbol`.
Step 156 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 157 [x]: Add or update targeted tests for `ValenceContext::new`, `ValenceContext::atom`, `ValenceContext::incident_bonds`, `periodic_table_row`, `required_valence_list`, and `rdkit_element_symbol`.
Step 158 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 159 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict valence`.

Step 160 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 161 [x]: Port `sanitize_cleanup_chirality_assignment` and `sanitize_adjust_hydrogens_assignment`.
Step 162 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 163 [x]: Add or update targeted tests for `sanitize_cleanup_chirality_assignment` and `sanitize_adjust_hydrogens_assignment`.
Step 164 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 165 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sanitize_adjust_hydrogens_assignment`.

Step 166 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 167 [x]: Audit the frozen baseline function list against the current code and rewrite `dev/gap_reports/smiles_rdkit_remaining_source_scan.md` with a zero-gap result if and only if every frozen function is fully covered.
Step 168 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 169 [x]: Update `dev/porting_inventory.md` and `crates/cosmolkit-core/src/support.rs` to reflect the final frozen-scope completion state.
Step 170 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 171 [x]: Run `cargo fmt --all`.
Step 172 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 173 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 174 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 175 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.
Step 176 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 177 [x]: Audit the frozen baseline against RDKit source one final time and write `dev/gap_reports/smiles_rdkit_frozen_coverage_final_audit.md` stating whether the frozen-scope coverage is exactly `100%`.

## Remaining Closure Plan

Step 178 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 179 [x]: Fix the currently failing strict tests by completing the corresponding RDKit-backed behaviors in `crates/cosmolkit-core/src/io/sdf.rs`, `crates/cosmolkit-core/src/notation/smiles.rs`, and `crates/cosmolkit-core/src/notation/smiles_write.rs` for `read_sdf_from_str_reads_simple_v3000_ctab`, `assign_chiral_types_from_bond_dirs_promotes_implicit_h_on_three_coordinate_center`, `handle_cx_part_and_name_clears_wiggly_single_bond_dirs_and_marks_unknown_stereo_like_rdkit`, `handle_cx_part_and_name_parses_double_bond_stereo_like_rdkit`, `mol_from_smiles_skip_cleanup_preserves_parser_temporary_props_like_rdkit_reader_flag`, `cx_individual_atom_labels_writes_label_entries`, `cx_individual_coords_writes_atom_order_coordinates`, `isomeric_writer_outputs_default_non_tetrahedral_chiral_classes_like_rdkit`, and `isomeric_writer_outputs_non_tetrahedral_chiral_classes_like_rdkit`.
Step 180 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 181 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict read_sdf_from_str_reads_simple_v3000_ctab assign_chiral_types_from_bond_dirs_promotes_implicit_h_on_three_coordinate_center handle_cx_part_and_name_clears_wiggly_single_bond_dirs_and_marks_unknown_stereo_like_rdkit handle_cx_part_and_name_parses_double_bond_stereo_like_rdkit mol_from_smiles_skip_cleanup_preserves_parser_temporary_props_like_rdkit_reader_flag cx_individual_atom_labels_writes_label_entries cx_individual_coords_writes_atom_order_coordinates isomeric_writer_outputs_default_non_tetrahedral_chiral_classes_like_rdkit isomeric_writer_outputs_non_tetrahedral_chiral_classes_like_rdkit`.

Step 182 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 183 [x]: Port and marker-close the remaining frozen parser functions in `crates/cosmolkit-core/src/notation/smiles.rs`: `assign_double_bond_stereo_from_directions`, `assign_chiral_types_from_bond_dirs`, `assign_chiral_types_from_3d`, `handle_cx_part_and_name`, `mol_from_smiles`, `to_mol`, `preprocess_smiles`, `parse_cx_extensions`, `parse_cx_doublebond_stereo`, `parse_cx_wedged_bonds`, `parse_cx_variable_attachments`, `parse_cx_coords`, `parse_cx_atom_labels`, `parse_cx_atom_values`, `parse_cx_atom_props`, `parse_cx_coordinate_bonds`, `parse_cx_zero_bonds`, `parse_cx_enhanced_stereo`, `parse_cx_unsaturation`, `parse_cx_ring_bonds`, `parse_cx_linknodes`, `parse_cx_data_sgroup_attr`, `parse_cx_data_sgroup`, `parse_cx_sgroup_hierarchy`, `parse_cx_polymer_sgroup`, `parse_cx_substitution`, `process_cx_radical_section`, `parse_cx_radicals`, and `remove_hs_update_explicit_count`, so that no `RDKit❌❌`, `RDKit❗❗`, `RDKit✔️❌`, `RDKit❗✔️`, or `RDKit✔️❗` marker remains inside the frozen parser scope.
Step 184 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 185 [x]: Add or update targeted tests covering the remaining parser closure work for the functions completed in Step 183.
Step 186 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 187 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict notation::smiles::tests`.

Step 188 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 189 [x]: Port and marker-close the remaining frozen writer functions in `crates/cosmolkit-core/src/notation/smiles_write.rs`: `mol_to_random_smiles_vect`, `mol_to_smiles_with_mode`, `prepare_plain_smiles_molecule`, `collect_fragment_write_plans`, `choose_fragment_start_atom`, `canonicalize_fragment_stack`, `write_mol_stack`, `mol_fragment_to_smiles`, `assemble_fragment_smiles`, `get_bond_smiles`, `get_molecule_bond_smiles`, `atom_needs_bracket`, `write_cx_zero_bonds`, `writer_parse_molfile_bond_cfg`, `write_cx_ringbond_cistrans_block`, `write_cx_linknodes_block`, `writer_parse_sgroup_crossings`, `cleanup_stereo_groups_for_cx_smiles`, `canonical_dfs_traversal`, `compute_writer_chiral_adjustments`, `insert_implicit_nontetrahedral_neighbors_for_writer`, `chiral_atom_needs_tag_inversion_for_writer`, `canonicalize_double_bond_for_writer`, `handle_dir_conflicts_across_double_bond_for_writer`, `remove_unwanted_bond_dir_specs_for_writer`, `remove_redundant_bond_dir_specs_for_writer`, and `same_side_dirs_are_compatible_for_writer`, so that no unresolved frozen-scope RDKit marker remains in the writer path.
Step 190 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 191 [x]: Add or update targeted tests covering the remaining writer closure work for the functions completed in Step 189.
Step 192 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 193 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict notation::smiles_write::tests`.

Step 194 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 195 [x]: Port and marker-close the remaining frozen canon-ranking functions in `crates/cosmolkit-core/src/notation/canon_rank.rs`: `rank_fragment_atoms`, `init_fragment_canon_atoms`, `init_canon_atoms_stereo_group_assignment`, `make_canon_bond_holder`, `rank_with_functor`, `refine_partitions`, `break_ties`, `hanoisort`, `hanoi`, `special_chirality_atom_compare`, `special_symmetry_atom_compare`, `atom_compare_base`, and `compare_ring_atoms_concerning_num_neighbors`, so that no unresolved frozen-scope RDKit marker remains in the canon-ranking path.
Step 196 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 197 [x]: Add or update targeted tests covering the remaining canon-ranking closure work for the functions completed in Step 195.
Step 198 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 199 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict canon_rank`.

Step 200 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 201 [x]: Port and marker-close the remaining frozen ring-perception functions in `crates/cosmolkit-core/src/chemistry/rings.rs`: `symmetrize_sssr_with_options`, `find_ring_families`, and `find_sssr_internal`, including the current non-URF and aggregate `findSSSR` or `symmetrizeSSSR` copied-source gaps identified by the audit, so that no unresolved frozen-scope RDKit marker remains in the ring path.
Step 202 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 203 [x]: Add or update targeted tests covering the remaining ring-perception closure work for the functions completed in Step 201.
Step 204 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 205 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict rings`.

Step 206 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 207 [x]: Port and marker-close the remaining frozen kekulization functions in `crates/cosmolkit-core/src/chemistry/kekulize.rs`: `kekulize_assignment`, `kekulize_fragment_assignment`, `kekulize_fused_assignment`, `mark_double_bond_candidates`, `kekulize_worker_matching`, `permute_dummies_and_match`, and `kekulize_if_possible_assignment`, so that no unresolved frozen-scope RDKit marker remains in the kekulize path.
Step 208 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 209 [x]: Add or update targeted tests covering the remaining kekulization closure work for the functions completed in Step 207.
Step 210 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 211 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize`.

Step 212 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 213 [x]: Port and marker-close the remaining frozen valence and sanitize functions in `crates/cosmolkit-core/src/chemistry/valence.rs` and `crates/cosmolkit-core/src/operations/ops.rs`: `ValenceContext::new`, `ValenceContext::atom`, `ValenceContext::incident_bonds`, `periodic_table_row`, `required_valence_list`, `rdkit_element_symbol`, `sanitize_cleanup_chirality_assignment`, and `sanitize_adjust_hydrogens_assignment`, and also close the remaining frozen sanitize/property/cleanup orchestration copied-source blocks in `operations/ops.rs`, so that no unresolved frozen-scope RDKit marker remains in the valence or sanitize path.
Step 214 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 215 [x]: Add or update targeted tests covering the remaining valence and sanitize closure work for the functions and orchestration completed in Step 213.
Step 216 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 217 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict valence sanitize`.

Step 218 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 219 [x]: Audit the frozen baseline against the current code and rewrite `dev/gap_reports/smiles_rdkit_remaining_source_scan.md` only if the entire frozen scope is now zero-gap; otherwise write the exact remaining function-level blockers and stop claiming closure.
Step 220 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 221 [x]: Update `dev/porting_inventory.md` and `crates/cosmolkit-core/src/support.rs` to reflect the new frozen-scope completion state without overstating parity or closure.
Step 222 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 223 [x]: Run `cargo fmt --all`.
Step 224 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 225 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 226 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 227 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.
Step 228 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 229 [x]: Audit the frozen baseline against RDKit source one final time and rewrite `dev/gap_reports/smiles_rdkit_frozen_coverage_final_audit.md`, and only state `Frozen-scope coverage is exactly 100%` if strict tests pass and no unresolved frozen-scope RDKit marker remains anywhere in the frozen baseline.

## Remaining Frozen Closure

Step 230 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 231 [x]: Read `dev/policy_invariants.md`, `dev/source_reproduction_protocol.md`, and `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`; record the exact frozen parser blockers that still apply before touching code.
Step 232 [x]: In `crates/cosmolkit-core/src/notation/smiles.rs`, audit `next_token` against the copied RDKit source lines and update only the markers that are proven by code and tests.
Step 233 [x]: Add or update focused tokenizer tests for the `next_token` branches audited in Step 232, including at least one negative or boundary token case when the source branch has one.
Step 234 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict notation::smiles::tests::token`.
Step 235 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 236 [x]: Audit and port `branch_open_token` in `smiles.rs`, including source comments, unsupported branches, marker status, and explicit error behavior.
Step 237 [x]: Add or update tests covering `branch_open_token`, unmatched branch handling, and branch error reporting.
Step 238 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict branch`.
Step 239 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 240 [x]: Audit and port `parse_mol` and `mol_from_smiles` together, because they share top-level parse orchestration and source error policy.
Step 241 [x]: Add or update tests covering top-level parse success, parse failure, source-unchanged behavior, and explicit unsupported behavior for Step 240.
Step 242 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol_from_smiles parse_mol`.
Step 243 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 244 [x]: Audit and port `finish_parse` and `cleanup_after_parsing`, including final sanitize dispatch, property cleanup, and error ordering against RDKit source.
Step 245 [x]: Add or update tests covering final parse cleanup, final sanitize dispatch, and cleanup error ordering for Step 244.
Step 246 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict finish_parse cleanup_after_parsing`.
Step 247 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 248 [x]: Audit and port `add_first_atom` and `add_disconnected_atom`, including disconnected fragment ordering, source atom properties, and marker accuracy.
Step 249 [x]: Add or update tests covering first-atom and disconnected-fragment parsing for Step 248.
Step 250 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict add_first_atom add_disconnected_atom`.
Step 251 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 252 [x]: Audit and port `add_atom_connected_to_active` and `add_branch_atom_connected_to_active`, including active atom updates and branch-local connection behavior.
Step 253 [x]: Add or update tests covering active atom and branch atom connection behavior for Step 252.
Step 254 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict add_atom_connected_to_active branch_atom`.
Step 255 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 256 [x]: Audit and port `add_single_bond_to_atom` and `add_branch_single_bond`, including implicit single-bond behavior and aromatic default handling.
Step 257 [x]: Add or update tests covering implicit single bonds, aromatic single-bond defaults, and branch single-bond behavior for Step 256.
Step 258 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict add_single_bond`.
Step 259 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 260 [x]: Audit and port `add_explicit_bond_to_atom` and `add_branch_explicit_bond`, including explicit bond type, direction, query-bond, and duplicate-bond error behavior.
Step 261 [x]: Add or update tests covering explicit atom bonds and explicit branch bonds for Step 260.
Step 262 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict explicit_bond`.
Step 263 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 264 [x]: Audit and port `add_ring_marker` and `add_single_bond_ring_marker`, including ring digit bookkeeping, duplicate ring marker errors, and implicit ring bond order behavior.
Step 265 [x]: Add or update tests covering single-bond ring markers and ring marker bookkeeping for Step 264.
Step 266 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ring_marker`.
Step 267 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 268 [x]: Audit and port `add_explicit_bond_ring_marker` and `close_ring_opening`, including explicit ring bond order, ring closure atom ordering, and error behavior.
Step 269 [x]: Add or update tests covering explicit ring closures, mismatched ring closure bonds, and ring closure error reporting for Step 268.
Step 270 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ring_closure`.
Step 271 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 272 [x]: Audit and port `check_ring_closure_branch_status`, `close_branch`, and `close_mol_rings`, including unclosed branch/ring diagnostics and final ring validation.
Step 273 [x]: Add or update tests covering unclosed branches, unclosed rings, branch closure order, and final ring closure diagnostics for Step 272.
Step 274 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict close_branch close_mol_rings`.
Step 275 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing parser work.
Step 276 [x]: Audit and port `add_frag_to_mol`, including fragment merge ordering, atom and bond id remapping, stereo group propagation, and source marker status.
Step 277 [x]: Add or update tests covering multi-fragment merge ordering and remapped atom or bond identities for Step 276.
Step 278 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict add_frag_to_mol fragment`.
Step 279 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing chirality parser work.
Step 280 [x]: Audit and port `adjust_atom_chirality_flags`, focusing only on tetrahedral atoms, ring-closure bond ordering, and tag inversion behavior.
Step 281 [x]: Add or update tests covering tetrahedral chirality adjustment, including the current fused-ring row-94 parse failure path.
Step 282 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict writer_canonical_fragment_scope_preserves_aromatic_fused_ring_form_like_rdkit_row_94`.
Step 283 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing chirality parser work.
Step 284 [x]: Audit and port `perturbation_order`, including RDKit `countSwapsToInterconvert` size and missing-element behavior without suppressing legitimate source errors.
Step 285 [x]: Add or update tests covering equal-size perturbation, size-mismatch failure, missing-probe-element failure, and the row-94 successful parse case if Step 280 depends on it.
Step 286 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict perturbation_order`.
Step 287 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing chirality parser work.
Step 288 [x]: Audit and port `chiral_atom_needs_tag_inversion`, `atom_has_fourth_valence`, and `is_unsaturated` together against RDKit source, including implicit-H/property-cache differences.
Step 289 [x]: Add or update tests covering explicit H, implicit-valence, query single-H, ring closure, and unsaturated-atom inversion branches for Step 288.
Step 290 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict chiral_atom_needs_tag_inversion atom_has_fourth_valence is_unsaturated`.
Step 291 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing stereo parser work.
Step 292 [x]: Audit and port `assign_stereochemistry_from_3d`, including conf-present and conf-null behavior, source copied lines, and explicit unsupported branches.
Step 293 [x]: Add or update tests covering 3D stereochemistry assignment with coordinates and the no-conformer path for Step 292.
Step 294 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict assign_stereochemistry_from_3d`.
Step 295 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing double-bond stereo parser work.
Step 296 [x]: Audit and port `set_double_bond_neighbor_directions`, including direction assignment from existing stereo, ring-direction handling, and source marker status.
Step 297 [x]: Add or update tests covering double-bond neighbor direction assignment from parsed stereo and existing stereo for Step 296.
Step 298 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict set_double_bond_neighbor_directions`.
Step 299 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing double-bond stereo parser work.
Step 300 [x]: Audit and port `set_bond_stereo_from_directions`, including cis/trans inference, unknown direction handling, and conflicting direction behavior.
Step 301 [x]: Add or update tests covering `BondStereo` inference from slash/backslash directions for Step 300.
Step 302 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict set_bond_stereo_from_directions`.
Step 303 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing query parser work.
Step 304 [x]: Audit and port `complete_smiles_query_scan_subset`, including every modeled query branch and explicit unsupported branches for unmodeled RDKit query behavior.
Step 305 [x]: Add or update tests covering the modeled query scan subset and unsupported query scan branches for Step 304.
Step 306 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict complete_smiles_query_scan_subset`.
Step 307 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before parser closure audit.
Step 308 [x]: Run a parser-only marker scan for `crates/cosmolkit-core/src/notation/smiles.rs` and list every remaining frozen-scope non-closed RDKit marker in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`.
Step 309 [x]: If Step 308 reports any parser blocker, add a new one-function-per-step continuation to this checklist before marking Step 309 complete; if it reports zero parser blockers, explicitly record zero parser blockers in the gap report.
Step 310 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict notation::smiles::tests`.
Step 311 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before kekulize closure work.
Step 312 [x]: Audit and port `kekulize_assignment` only, including source markers, explicit unsupported behavior, and parity fixtures.
Step 313 [x]: Add or update tests covering `kekulize_assignment` success and failure branches for Step 312.
Step 314 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_assignment`.
Step 315 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing kekulize closure work.
Step 316 [x]: Audit and port `kekulize_fragment_assignment` only, including fragment bitset behavior, source markers, and error propagation.
Step 317 [x]: Add or update tests covering `kekulize_fragment_assignment` branch coverage for Step 316.
Step 318 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_fragment_assignment`.
Step 319 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing kekulize closure work.
Step 320 [x]: Audit and port `kekulize_fused_assignment` only, including fused ring grouping, candidate propagation, and marker status.
Step 321 [x]: Add or update tests covering fused-ring kekulization branches for Step 320.
Step 322 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_fused_assignment fused`.
Step 323 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing kekulize closure work.
Step 324 [x]: Audit and port `mark_double_bond_candidates` only, including every candidate exclusion branch.
Step 325 [x]: Add or update tests covering double-bond candidate inclusion and exclusion branches for Step 324.
Step 326 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mark_double_bond_candidates`.
Step 327 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing kekulize closure work.
Step 328 [x]: Audit and port `kekulize_worker_matching` only, including matching success, no-match, timeout, and error branches.
Step 329 [x]: Add or update tests covering `kekulize_worker_matching` branches for Step 328.
Step 330 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_worker_matching`.
Step 331 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing kekulize closure work.
Step 332 [x]: Audit and port `permute_dummies_and_match` only, including dummy permutation order and failure handling.
Step 333 [x]: Add or update tests covering dummy permutation and matching branches for Step 332.
Step 334 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict permute_dummies_and_match`.
Step 335 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before completing kekulize closure work.
Step 336 [x]: Audit and port `kekulize_if_possible_assignment` only, including success, fallback, and explicit failure branches.
Step 337 [x]: Add or update tests covering `kekulize_if_possible_assignment` branches for Step 336.
Step 338 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_if_possible_assignment`.
Step 339 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before kekulize closure audit.
Step 340 [x]: Run a kekulize-only marker scan for `crates/cosmolkit-core/src/chemistry/kekulize.rs` and list every remaining frozen-scope non-closed RDKit marker in `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`.
Step 341 [x]: If Step 340 reports any kekulize blocker, add a new one-function-per-step continuation to this checklist before marking Step 341 complete; if it reports zero kekulize blockers, explicitly record zero kekulize blockers in the gap report.
Step 342 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize`.
Step 343 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before valence closure work.
Step 344 [x]: Audit and port `ValenceContext::new`, `ValenceContext::atom`, and `ValenceContext::incident_bonds` only, including boundary and invalid-id behavior.
Step 345 [x]: Add or update tests covering `ValenceContext` construction and access branches for Step 344.
Step 346 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ValenceContext`.
Step 347 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing valence closure work.
Step 348 [x]: Audit and port `periodic_table_row`, `required_valence_list`, and `rdkit_element_symbol` only, including unknown element and edge periodic-table branches.
Step 349 [x]: Add or update tests covering element row, valence list, and symbol branches for Step 348.
Step 350 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict periodic_table_row required_valence_list rdkit_element_symbol`.
Step 351 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before sanitize closure work.
Step 352 [x]: Audit and port `sanitize_cleanup_chirality_assignment` only, including tetrahedral, non-tetrahedral, stereo group, and invalid permutation branches.
Step 353 [x]: Add or update tests covering cleanup chirality branches for Step 352.
Step 354 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sanitize_cleanup_chirality_assignment`.
Step 355 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing sanitize closure work.
Step 356 [x]: Audit and port `sanitize_adjust_hydrogens_assignment` only, including pyrrolic H, explicit H, valence interaction, and no-op branches.
Step 357 [x]: Add or update tests covering adjust-hydrogens branches for Step 356.
Step 358 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sanitize_adjust_hydrogens_assignment`.
Step 359 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before sanitize orchestration closure work.
Step 360 [x]: Audit remaining sanitize/property/cleanup orchestration copied-source blocks in `crates/cosmolkit-core/src/operations/ops.rs`, one block at a time, and update only markers proven by implementation and tests.
Step 361 [x]: Add or update tests for every orchestration branch changed in Step 360.
Step 362 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict operations::ops::tests::sanitized`.
Step 363 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before valence/sanitize closure audit.
Step 364 [x]: Run a valence/sanitize marker scan for `crates/cosmolkit-core/src/chemistry/valence.rs` and `crates/cosmolkit-core/src/operations/ops.rs`; write exact remaining blockers to `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`.
Step 365 [x]: If Step 364 reports any valence or sanitize blocker, add a new one-function-or-one-block-per-step continuation to this checklist before marking Step 365 complete; if it reports zero blockers, explicitly record zero valence/sanitize blockers in the gap report.
Step 366 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict valence sanitize`.
Step 367 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before final marker audit.
Step 368 [x]: Run a full frozen-scope RDKit marker scan across parser, writer, canon ranking, rings, kekulize, valence, sanitize, and operation orchestration files; write the complete result to `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`.
Step 369 [x]: If Step 368 reports any blocker, append a new fine-grained continuation section to this checklist before marking Step 369 complete; do not claim frozen-scope closure.
Step 370 [ ]: If Step 368 reports zero blockers, update `dev/porting_inventory.md` and `crates/cosmolkit-core/src/support.rs` to reflect the exact frozen-scope state without overstating parity.
Step 371 [x]: Run `cargo fmt --all`.
Step 372 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before strict validation.
Step 373 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 374 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before final strict tests.
Step 375 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict notation::smiles::tests notation::smiles_write::tests canon_rank rings kekulize valence sanitize`.
Step 376 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before full strict validation.
Step 377 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.
Step 378 [ ]: If Step 377 fails, write each failing test, root cause, and required one-step continuation to `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`; do not claim final closure.
Step 379 [ ]: If Step 377 passes and Step 368 reports zero blockers, audit the frozen baseline against RDKit source one final time and rewrite `dev/gap_reports/smiles_rdkit_frozen_coverage_final_audit.md`.
Step 380 [x]: Only after Steps 368, 377, and 379 all pass, state `Frozen-scope coverage is exactly 100%`; otherwise state the exact remaining blockers and leave the checklist open.

## Kekulize Continuation After Step 340 Scan

Step 381 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing the remaining `kekulize_fragment_assignment` closure blocks.
Step 382 [x]: Audit and port only the ring-info, wedged-ring ordering, and filtered ring-copy copied-source block inside `kekulize_fragment_assignment`.
Step 383 [x]: Add or update focused tests covering the `kekulize_fragment_assignment` ring-ordering and filtered-ring branches for Step 382.
Step 384 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict ordered_kekulize_rings kekulize_fragment`.
Step 385 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing the remaining `kekulize_fragment_assignment` closure blocks.
Step 386 [x]: Audit and port only the aromatic-flag clearing and pyrrolic explicit-H reset copied-source block inside `kekulize_fragment_assignment`.
Step 387 [x]: Add or update focused tests covering the `kekulize_fragment_assignment` aromatic-flag clearing and explicit-H reset branches for Step 386.
Step 388 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_fragment clear_is_disabled preserves_aromatic_flags`.
Step 389 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing kekulize worker closure work.
Step 390 [x]: Audit and port `kekulize_worker_matching` copied-source block to closure-level evidence, updating only markers proven by code shape and tests.
Step 391 [x]: Add or update focused tests covering the remaining `kekulize_worker_matching` option ordering, wedge-priority, rollback, and terminal failure branches for Step 390.
Step 392 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict kekulize_worker_matching worker_sorted_atoms back_track`.
Step 393 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing dummy-permutation closure work.
Step 394 [x]: Audit and port `permute_dummies_and_match` copied-source blocks, including `QuestionEnumerator`, permutation order, and reset behavior.
Step 395 [x]: Add or update focused tests covering `permute_dummies_and_match` mask order, no-question exit, no-backtrack recovery, and final failure handling for Step 394.
Step 396 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict permute_dummies_and_match question_switch_masks`.

## Valence/Sanitize Continuation After Step 364 Scan

Step 397 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `ValenceContext::new` closure work.
Step 398 [x]: Audit and port only the remaining performance-gap block in `ValenceContext::new`.
Step 399 [x]: Add or update focused tests covering the remaining `ValenceContext::new` cache-present and cache-absent closure evidence for Step 398.
Step 400 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict valence_context_`.

Step 401 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `add_hs_terminal_coords` closure work.
Step 402 [x]: Audit and port only the remaining copied-source block in `add_hs_terminal_coords`.
Step 403 [x]: Add or update focused tests covering `add_hs_terminal_coords` degree and geometry branches for Step 402.
Step 404 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict add_hs_terminal_coords`.

Step 405 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_after_remove_hs_removal` closure work.
Step 406 [x]: Audit and port only the remaining copied-source block in `sanitize_after_remove_hs_removal`.
Step 407 [x]: Add or update focused tests covering `sanitize_after_remove_hs_removal` sanitize-trigger branches for Step 406.
Step 408 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict sanitize_after_remove_hs_removal remove_hs`.

Step 409 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `run_sanitize_pipeline` closure work.
Step 410 [x]: Audit and port only the remaining sanitize orchestration copied-source blocks in `run_sanitize_pipeline`.
Step 411 [x]: Add or update focused tests covering the remaining `run_sanitize_pipeline` stage-order and cache-update branches for Step 410.
Step 412 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict operations::ops::tests::sanitized`.

Step 413 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_nitrogens_cleanup_assignment` closure work.
Step 414 [x]: Audit and port only the remaining copied-source block in `sanitize_nitrogens_cleanup_assignment`.
Step 415 [x]: Add or update focused tests covering the remaining nitro cleanup selection and rewrite branches for Step 414.
Step 416 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict neutral_nitro nitrogens_cleanup`.

Step 417 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_phosphorus_cleanup_assignment` closure work.
Step 418 [x]: Audit and port only the remaining copied-source block in `sanitize_phosphorus_cleanup_assignment`.
Step 419 [x]: Add or update focused tests covering the remaining phosphorus oxo cleanup branches for Step 418.
Step 420 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict phosphorus_oxo phosphorus_cleanup`.

Step 421 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_halogen_cleanup_assignment` closure work.
Step 422 [x]: Audit and port only the remaining copied-source block in `sanitize_halogen_cleanup_assignment`.
Step 423 [x]: Add or update focused tests covering the remaining hypervalent halogen cleanup branches for Step 422.
Step 424 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict halogen_oxo halogen_cleanup`.

Step 425 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_cleanup_incident_bonds` closure work.
Step 426 [x]: Audit and port only the remaining copied-source block in `sanitize_cleanup_incident_bonds`.
Step 427 [x]: Add or update focused tests covering `sanitize_cleanup_incident_bonds` adjacency and bond-order branches for Step 426.
Step 428 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict cleanup_incident_bonds`.

Step 429 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_organometallic_cleanup_assignment` closure work.
Step 430 [x]: Audit and port only the remaining copied-source block in `sanitize_organometallic_cleanup_assignment`.
Step 431 [x]: Add or update focused tests covering the remaining organometallic cleanup ordering and early-exit branches for Step 430.
Step 432 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict organometallic_cleanup`.

Step 433 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_metal_bond_cleanup_assignment` closure work.
Step 434 [x]: Audit and port only the remaining copied-source block in `sanitize_metal_bond_cleanup_assignment`.
Step 435 [x]: Add or update focused tests covering the remaining metal-bond cleanup ranking and bond-rewrite branches for Step 434.
Step 436 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict metal_bond_cleanup`.

Step 437 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_is_hypervalent_nonmetal` closure work.
Step 438 [x]: Audit and port only the remaining copied-source block in `sanitize_is_hypervalent_nonmetal`.
Step 439 [x]: Add or update focused tests covering the remaining hypervalent-nonmetal predicate branches for Step 438.
Step 440 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict hypervalent_nonmetal`.

Step 441 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before continuing `sanitize_organometallic_single_bonded_metals` closure work.
Step 442 [x]: Audit and port only the remaining copied-source block in `sanitize_organometallic_single_bonded_metals`.
Step 443 [x]: Add or update focused tests covering the remaining single-bonded metal selection branches for Step 442.
Step 444 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict single_bonded_metals organometallic`.

## Frozen-Scope Continuation After Step 368 Scan

Step 445 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before resuming frozen-scope blocker closure after the Step 368 full scan.
Step 446 [x]: Complete the remaining kekulize continuation block in Steps 381-396 before claiming any new frozen-scope closure result.
Step 447 [x]: Complete the remaining valence/sanitize continuation block in Steps 397-444 before claiming any new frozen-scope closure result.
Step 448 [x]: After Steps 381-444 finish, rerun the full frozen-scope marker scan and append any further one-function-or-one-block continuation required by the new result before claiming frozen-scope closure.

## Frozen-Scope Continuation After Step 448 Scan

Step 449 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before closing the remaining frozen-scope blockers revealed by the Step 448 full scan.
Step 450 [x]: Audit and port only the remaining performance-gap block in `ValenceContext::new` by preserving immutable topology adjacency on constructed molecules.
Step 451 [x]: Add or update focused tests covering topology-adjacency reuse and empty-topology normalization for Step 450.
Step 452 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict valence_context_`.
Step 453 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` before closing the remaining `add_hs_terminal_coords` copied-source block.
Step 454 [x]: Audit and port only the remaining copied-source block in `add_hs_terminal_coords`.
Step 455 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict add_hs_terminal_coords`.
Step 456 [x]: Rerun the full frozen-scope marker scan plus strict validation (`cargo check -p cosmolkit-core --features op-contracts-strict` and `cargo test -p cosmolkit-core --features op-contracts-strict`) and record the zero-blocker result.
