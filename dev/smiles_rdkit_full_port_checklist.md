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

Step 212 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 213 [ ]: Port and marker-close the remaining frozen valence and sanitize functions in `crates/cosmolkit-core/src/chemistry/valence.rs` and `crates/cosmolkit-core/src/operations/ops.rs`: `ValenceContext::new`, `ValenceContext::atom`, `ValenceContext::incident_bonds`, `periodic_table_row`, `required_valence_list`, `rdkit_element_symbol`, `sanitize_cleanup_chirality_assignment`, and `sanitize_adjust_hydrogens_assignment`, and also close the remaining frozen sanitize/property/cleanup orchestration copied-source blocks in `operations/ops.rs`, so that no unresolved frozen-scope RDKit marker remains in the valence or sanitize path.
Step 214 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 215 [ ]: Add or update targeted tests covering the remaining valence and sanitize closure work for the functions and orchestration completed in Step 213.
Step 216 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 217 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict valence sanitize`.

Step 218 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 219 [ ]: Audit the frozen baseline against the current code and rewrite `dev/gap_reports/smiles_rdkit_remaining_source_scan.md` only if the entire frozen scope is now zero-gap; otherwise write the exact remaining function-level blockers and stop claiming closure.
Step 220 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 221 [ ]: Update `dev/porting_inventory.md` and `crates/cosmolkit-core/src/support.rs` to reflect the new frozen-scope completion state without overstating parity or closure.
Step 222 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 223 [ ]: Run `cargo fmt --all`.
Step 224 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 225 [ ]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 226 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 227 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.
Step 228 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 229 [ ]: Audit the frozen baseline against RDKit source one final time and rewrite `dev/gap_reports/smiles_rdkit_frozen_coverage_final_audit.md`, and only state `Frozen-scope coverage is exactly 100%` if strict tests pass and no unresolved frozen-scope RDKit marker remains anywhere in the frozen baseline.
