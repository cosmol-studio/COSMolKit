# RDKit MOL2 Full Port Plan

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
- Copying C++ comments, adding a dispatch stub, or adding placeholder branches is not a completed `Port` step.
- Do not use "smallest subpart", skeleton code, dispatch-only code, placeholder code, TODO-only code, or partial porting as completion.
- If a step is too large to complete fully, regenerate the plan by splitting that step into smaller full-port steps.
- Each split step must still represent a complete source-backed behavior, not a placeholder.
- A `Port` step is complete only when the selected RDKit behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Source Authority

This plan targets RDKit-compatible MOL2 parsing, not an independent Tripos-format parser.

Primary source anchors:

- `third_party/rdkit/Code/GraphMol/FileParsers/FileParsers.h`
- `third_party/rdkit/Code/GraphMol/FileParsers/Mol2FileParser.cpp`
- `third_party/rdkit/Code/GraphMol/FileParsers/testMol2ToMol.cpp`
- `third_party/rdkit/Code/GraphMol/FileParsers/test_data/*.mol2`

The local repository currently contains RDKit MOL2 source and tests but no OpenBabel source tree. Do not use OpenBabel behavior unless its exact source is vendored or otherwise approved as a source authority.

## Target Rust Artifacts

- `dev/gap_reports/rdkit_mol2_source_map.md`
- `crates/cosmolkit-core/src/io/mol2.rs`
- `crates/cosmolkit-core/src/io/mod.rs`
- `crates/cosmolkit-core/src/support.rs`
- `crates/cosmolkit-core/src/lib.rs`
- `crates/cosmolkit-core/tests/rdkit_mol2_read_parity.rs`
- `tests/golden/mol2_read.jsonl`
- `tests/fixtures/mol2/`

## Target Rust Function Surface

- `Mol2ReadError`
- `Mol2Type`
- `Mol2ReadParams`
- `Mol2Record`
- `read_mol2_from_str`
- `read_mol2_from_str_with_params`
- `read_mol2_file`
- `read_mol2_file_with_params`
- `mol_from_mol2_data_stream_like_rdkit`
- `mol_from_mol2_block_like_rdkit`
- `mol_from_mol2_file_like_rdkit`
- `scan_mol2_sections_like_rdkit`
- `parse_mol2_molecule_header_like_rdkit`
- `parse_mol2_atom_line_like_rdkit`
- `parse_mol2_bond_line_like_rdkit`
- `parse_mol2_atom_block_like_rdkit`
- `parse_mol2_bond_block_like_rdkit`
- `read_formal_charges_from_attr_like_rdkit`
- `guess_formal_charges_like_rdkit`
- `fix_nitro_substructure_and_charge_like_rdkit`
- `check_no_h_neighbors_n_oxide_like_rdkit`
- `clean_up_mol2_substructures_like_rdkit`
- `assign_chiral_types_from_3d_for_mol2_like_rdkit`
- `sanitize_mol2_molecule_like_rdkit`
- `finish_mol2_read_like_rdkit`

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit current COSMolKit MOL2 support and write `dev/gap_reports/rdkit_mol2_source_map.md` with current-state evidence.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Audit RDKit `FileParsers.h` MOL2 declarations and append the Rust API mapping to `dev/gap_reports/rdkit_mol2_source_map.md`.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Audit RDKit `Mol2FileParser.cpp` function bodies and append the function-by-function source map to `dev/gap_reports/rdkit_mol2_source_map.md`.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Audit RDKit `testMol2ToMol.cpp` and append the parity fixture matrix to `dev/gap_reports/rdkit_mol2_source_map.md`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Update `crates/cosmolkit-core/src/support.rs` with `MOL2_READ_FEATURE` marked unsupported until source-backed implementation steps complete.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Implement `Mol2ReadError`, `Mol2Type`, `Mol2ReadParams`, and `Mol2Record` in `crates/cosmolkit-core/src/io/mol2.rs` from RDKit `Mol2Type` and `Mol2ParserParams`.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Port the section-scanning portion of RDKit `MolFromMol2DataStream` as `scan_mol2_sections_like_rdkit`.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Add MOL2 section-scanning tests for missing `MOLECULE`, missing `ATOM`, repeated `MOLECULE`, `BOND`, and `UNITY_ATOM_ATTR` boundaries.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_section -- --nocapture`.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Port the molecule-name, counts-line, molecule-type skip, and charge-type parsing portion of RDKit `MolFromMol2DataStream` as `parse_mol2_molecule_header_like_rdkit`.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Add MOL2 molecule-header tests for name trimming, empty counts line, invalid atom count, missing bond count default, zero atoms, molecule type ignore, and charge-type property retention.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_header -- --nocapture`.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Port RDKit `ParseMol2FileAtomLine` atom-id skip, atom-name property, coordinate parsing, SYBYL type property, no-implicit-hydrogen policy, and partial-charge property as `parse_mol2_atom_line_like_rdkit`.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Add MOL2 atom-line tests for required fields, coordinate errors, `_TriposAtomName`, `_TriposAtomType`, `_TriposPartialCharge`, and no-implicit-hydrogen state.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_atom_line -- --nocapture`.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port RDKit `ParseMol2FileAtomLine` LP, ANY, Du, HEV, HET, HAL, and periodic-table symbol handling into `parse_mol2_atom_line_like_rdkit`.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Add MOL2 atom-type tests for LP removal, ANY query, Du query, HEV query, HET query, HAL query, ordinary symbols, and unknown-symbol failure.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_atom_type -- --nocapture`.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Port RDKit `ParseMol2AtomBlock` as `parse_mol2_atom_block_like_rdkit`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Add MOL2 atom-block tests for LP index correspondence, 3D conformer creation, explicit-hydrogen detection, premature EOF, and atom-count postcondition.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_atom_block -- --nocapture`.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Port RDKit `ParseMol2FileBondLine` as `parse_mol2_bond_line_like_rdkit`.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Add MOL2 bond-line tests for one-based indices, LP-skipped endpoints, index mismatch, bond types `1`, `am`, `2`, `3`, `ar`, `du`, `un`, `nc`, and unsupported strings.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_bond_line -- --nocapture`.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Port RDKit `ParseMol2BondBlock` as `parse_mol2_bond_block_like_rdkit`.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Add MOL2 bond-block tests for aromatic bond flags, aromatic atom flags, skipped bad bonds, premature EOF, and bond-count postcondition.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_bond_block -- --nocapture`.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Port RDKit `readFormalChargesFromAttr` as `read_formal_charges_from_attr_like_rdkit`.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Add MOL2 `UNITY_ATOM_ATTR` tests for `AtomExpr` charge assignment, malformed header, malformed charge, section termination, blank-line termination, comment termination, and premature EOF.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_unity_atom_attr -- --nocapture`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Port RDKit `fixNitroSubstructureAndCharge` as `fix_nitro_substructure_and_charge_like_rdkit`.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Add MOL2 nitro-cleanup tests for two double-bond oxygen neighbors and non-matching nitrogen environments.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_nitro_cleanup -- --nocapture`.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Port RDKit `chkNoHNeighbNOx` as `check_no_h_neighbors_n_oxide_like_rdkit`.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Add MOL2 N-oxide helper tests for hydrogen-neighbor counting and terminal oxygen detection.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_n_oxide_helper -- --nocapture`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Port RDKit `cleanUpMol2Substructures` `N.4` and `O.co2` branches as `clean_up_mol2_substructures_like_rdkit`.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Add MOL2 cleanup tests for `N.4`, carboxylate `O.co2`, sulfonate `O.co2`, phosphate `O.co2`, degree errors, and unsupported neighbor errors.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_cleanup_o_co2 -- --nocapture`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Port RDKit `cleanUpMol2Substructures` two-nitrogen `C.cat` branch as `clean_up_mol2_substructures_like_rdkit`.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Add MOL2 cleanup tests for two-nitrogen `C.cat` N-oxide precedence, hydrogen-count precedence, ring precedence, and fallback precedence.
Step 84 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 85 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_cleanup_c_cat_two_n -- --nocapture`.
Step 86 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 87 [x]: Port RDKit `cleanUpMol2Substructures` three-nitrogen `C.cat` branch as `clean_up_mol2_substructures_like_rdkit`.
Step 88 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 89 [x]: Add MOL2 cleanup tests for three-nitrogen `C.cat`, bad nitrogen counts, already-fixed neighbors, and heavy-atom-degree tie handling.
Step 90 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 91 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_cleanup_c_cat_three_n -- --nocapture`.
Step 92 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 93 [x]: Port RDKit `guessFormalCharges` as `guess_formal_charges_like_rdkit`.
Step 94 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 95 [x]: Add MOL2 formal-charge guessing tests for non-carbon atoms, query atoms, aromatic five-membered rings, `N.ar` three-aromatic-bond warning behavior, multi-valence atoms, aromatic charge clamping, and nitro repair handoff.
Step 96 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 97 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_guess_charges -- --nocapture`.
Step 98 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 99 [x]: Port the unsanitized assembly path of RDKit `MolFromMol2DataStream` as `mol_from_mol2_data_stream_like_rdkit`.
Step 100 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 101 [x]: Add MOL2 unsanitized assembly tests for atom block, bond block, no bonds, cleanup enabled, cleanup disabled, `UNITY_ATOM_ATTR` charge override, and null cleanup result.
Step 102 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 103 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_data_stream_unsanitized -- --nocapture`.
Step 104 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 105 [x]: Port RDKit `MolOps::assignChiralTypesFrom3D` invocation semantics as `assign_chiral_types_from_3d_for_mol2_like_rdkit`.
Step 106 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 107 [x]: Add MOL2 3D chirality tests for first-conformer assignment before hydrogen removal.
Step 108 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 109 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_chirality_3d -- --nocapture`.
Step 110 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 111 [x]: Port RDKit sanitized-with-removeHs branch as `sanitize_mol2_molecule_like_rdkit`.
Step 112 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 113 [x]: Add MOL2 sanitize tests for cleanup before stereo detection, bond stereochemistry before hydrogen removal, hydrogen removal without resanitize, skipped organometallic cleanup, final sanitize, property-cache update, and stereochemistry assignment.
Step 114 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 115 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_sanitize_remove_hs -- --nocapture`.
Step 116 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 117 [x]: Port RDKit sanitized-without-removeHs branch as `sanitize_mol2_molecule_like_rdkit`.
Step 118 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 119 [x]: Add MOL2 sanitize tests for no hydrogen removal, skipped organometallic cleanup, bond-stereo detection after sanitize, property-cache update, and stereochemistry assignment.
Step 120 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 121 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_sanitize_keep_hs -- --nocapture`.
Step 122 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 123 [x]: Implement `finish_mol2_read_like_rdkit` to apply the exact RDKit ordering from `MolFromMol2DataStream`.
Step 124 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 125 [ ]: Add end-to-end MOL2 unit tests for `pyrazole_pyridine.mol2`, `benzene.mol2`, `mol_noatoms.mol2`, `mol_nomol.mol2`, `lonePairMol.mol2`, `symmetricGuanidine.mol2`, `Noxide.mol2`, `fusedRing.mol2`, and `EZ_mol2_issue114.mol2`.
Step 126 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 127 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_end_to_end -- --nocapture`.
Step 128 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 129 [x]: Port `Mol2BlockToMol`, `Mol2FileToMol`, `MolFromMol2Block`, and `MolFromMol2File` wrappers as `read_mol2_from_str`, `read_mol2_from_str_with_params`, `read_mol2_file`, `read_mol2_file_with_params`, `mol_from_mol2_block_like_rdkit`, and `mol_from_mol2_file_like_rdkit`.
Step 130 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 131 [ ]: Add wrapper-level MOL2 unit tests for default params, file errors, empty-file behavior, and block dispatch.
Step 132 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 133 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict mol2_wrapper -- --nocapture`.
Step 134 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 135 [x]: Add RDKit-generated `tests/golden/mol2_read.jsonl` covering topology, atom fields, bond fields, charges, aromaticity, 3D coordinates, chirality, and canonical SMILES for the MOL2 fixture matrix.
Step 136 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 137 [x]: Run `cargo test -p cosmolkit-core --features op-contracts-strict --test rdkit_mol2_read_parity -- --nocapture`.
Step 138 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 139 [x]: Modify `crates/cosmolkit-core/src/io/mod.rs` and `crates/cosmolkit-core/src/lib.rs` to expose the completed MOL2 reader surface.
Step 140 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 141 [ ]: Update `crates/cosmolkit-core/src/support.rs` to mark `MOL2_READ_FEATURE` experimental with an exact support description and remaining unsupported behavior if any remains.
Step 142 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 143 [ ]: Update `dev/porting_inventory.md` and `tests/README.md` with the MOL2 reader status, fixture notes, and parity boundary.
Step 144 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 145 [ ]: Run `cargo fmt --all`.
Step 146 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 147 [ ]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 148 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 149 [ ]: Run `cargo test -p cosmolkit-core --features op-contracts-strict`.
Step 150 [ ]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 151 [ ]: Audit `crates/cosmolkit-core/src/io/mol2.rs` markers and update `dev/gap_reports/rdkit_mol2_source_map.md` with final marker counts and any explicit unsupported branches.
