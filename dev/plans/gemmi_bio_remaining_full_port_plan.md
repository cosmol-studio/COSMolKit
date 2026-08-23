# Gemmi BioStructure Remaining Full Port Plan

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
- A `Port` step is complete only when the selected Gemmi behavior is fully implemented, source-backed, tested, and no placeholder remains.
- During porting, do not execute any git operation unless the user explicitly asks for it.
- This includes read-only git operations such as `git status`, `git diff`, and `git log`.

## Source Baseline

- Gemmi commit: `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`.
- Primary writer source: `third_party/gemmi/src/to_mmcif.cpp`.
- Output-group source: `third_party/gemmi/include/gemmi/to_mmcif.hpp`.
- CIF serializer source: `third_party/gemmi/include/gemmi/to_cif.hpp`.
- CIF document mutation and quoting source: `third_party/gemmi/include/gemmi/cifdoc.hpp`.
- Public ownership boundary: structural writers belong to `BioStructure`, not `Protein` or `Molecule`.

## Plan

Step 1 [x]: Read `dev/agent_plan_standard.md`.
Step 2 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 3 [x]: Audit every remaining BioStructure implementation gap against the pinned Gemmi source and write the function-and-field-level report to `dev/gap_reports/gemmi_bio_remaining_full_port_inventory.md`.
Step 4 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 5 [x]: Port the Gemmi CIF value quoting, item-category replacement, loop construction, and block mutation primitives required by `update_mmcif_block` into the single internal CIF document implementation.
Step 6 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 7 [x]: Add focused tests for CIF value quoting, category replacement, pair/loop replacement, row-width enforcement, and stable item ordering.
Step 8 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 9 [x]: Run `cargo test -p cosmolkit-core --lib io::bio::tests::cif_writer_ --features op-contracts-strict`.
Step 10 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 11 [x]: Port Gemmi `cif::WriteOptions`, text-field writing, pair writing, loop writing, category separation, block writing, and document writing into the internal CIF serializer.
Step 12 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 13 [x]: Add focused exact-output tests for every CIF serialization option and multiline or CRLF text-field boundary.
Step 14 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 15 [x]: Run `cargo test -p cosmolkit-core --lib io::bio::tests::cif_serializer_ --features op-contracts-strict`.
Step 16 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 17 [x]: Port `MmcifOutputGroups`, numeric/null formatting helpers, block-name validation, and `write_cell_parameters` with a public options type that preserves Gemmi defaults.
Step 18 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 19 [x]: Port `add_cif_atoms` including optional columns, group-PDB policy, auth columns, model numbering, anisotropic rows, calc flags, TLS identifiers, deuterium fractions, and allocation behavior.
Step 20 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 21 [x]: Add focused tests for atom-site and anisotropic writer categories across single-model, multi-model, optional-column, alternate-location, charge, and null-value branches.
Step 22 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 23 [x]: Run `cargo test -p cosmolkit-core --lib io::bio::writer::tests::mmcif_atom_writer_ --features op-contracts-strict`.
Step 24 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 25 [x]: Port `xmeric_to_number` and `write_assemblies` including operator deduplication, transform serialization, subchain expansion, properties, and oligomeric-count behavior.
Step 26 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 27 [x]: Port `write_ncs_oper`, `write_struct_conn`, and `write_cispeps` including identity insertion, partner lookup, symmetry and distance fields, connection-type categories, model selection, and missing-partner omission.
Step 28 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 29 [x]: Add focused tests for assemblies, NCS operations, structural connections, and cis-peptide categories including omission and deduplication branches.
Step 30 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 31 [x]: Run `cargo test -p cosmolkit-core --lib io::bio::writer::tests::mmcif_structural_category_writer_ --features op-contracts-strict`.
Step 32 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 33 [x]: Port the entry, database-status, author, cell, symmetry, entity, entity-polymer, structure-reference, chemical-component, experiment, diffraction, reflection, refinement, title, and keyword branches of `update_mmcif_block` for every state represented by `BioStructure`.
Step 34 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 35 [x]: Port the struct-asym, ORIGX, helix, sheet, biological-detail, modified-residue, SCALE, atom-type, entity-sequence, TLS, and software branches of `update_mmcif_block` for every state represented by `BioStructure`.
Step 36 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 37 [x]: Add focused tests for every `MmcifOutputGroups` branch and every represented metadata category, including category suppression and replacement in a pre-populated block.
Step 38 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 39 [x]: Run `cargo test -p cosmolkit-core --lib io::bio::tests::mmcif_output_group_ --features op-contracts-strict`.
Step 40 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 41 [x]: Port `update_mmcif_block`, `make_mmcif_document`, `make_mmcif_block`, `make_mmcif_headers`, and `add_minimal_mmcif_data` as the single internal structural-writer call path.
Step 42 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 43 [x]: Add `BioStructure::to_mmcif`, `BioStructure::to_mmcif_with_options`, and `BioStructure::write_mmcif` only after the complete writer path exists, with typed I/O errors and no Protein or Molecule aliases.
Step 44 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 45 [x]: Add public Rust integration tests for PDB-to-mmCIF writing, mmCIF structural roundtrip, file writing, option groups, hierarchy preservation, coordinates, identifiers, metadata, assemblies, connections, and multi-model structures.
Step 46 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 47 [x]: Run `cargo test -p cosmolkit-core --test bio_mmcif_writer --release --features op-contracts-strict`.
Step 48 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 49 [x]: Add reproducible pinned-Gemmi golden-data generation tooling and committed small-fixture provenance for exact category, ordering, quoting, and serialization comparison.
Step 50 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 51 [x]: Run the Gemmi writer golden-data preparation command documented by Step 49.
Step 52 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 53 [x]: Add public exact-output and semantic-reparse parity tests against every prepared pinned-Gemmi writer fixture and lock all discovered edge cases locally.
Step 54 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 55 [x]: Run `cargo test -p cosmolkit-core --test bio_mmcif_gemmi_parity --release --features op-contracts-strict`.
Step 56 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 57 [x]: Expose the completed BioStructure mmCIF writer and options through Python with Rust-equivalent ownership, errors, defaults, and no Protein or Molecule writer aliases.
Step 58 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 59 [x]: Regenerate Python stubs and add Python tests for string writing, file writing, options, exception conversion, source immutability, and structural reparse.
Step 60 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 61 [x]: Run `.venv/bin/pytest python/tests -q -k 'biostructure and mmcif'`.
Step 62 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 63 [x]: Close the remaining non-green Gemmi reader markers and strict-mode guards when their source behavior is representable, and document any genuinely unmodeled state without weakening an existing supported boundary.
Step 64 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 65 [x]: Add focused regression tests for every reader marker or strict-mode branch changed in Step 63.
Step 66 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 67 [x]: Run `cargo test -p cosmolkit-core --lib io::bio::tests --release --features op-contracts-strict`.
Step 68 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 69 [x]: Update current BioStructure I/O policy, public Rust/Python documentation, examples, support matrices, and the wwPDB stress protocol to describe the implemented writer boundary without editing historical plans or archives.
Step 70 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 71 [x]: Audit the completed Bio line against the pinned Gemmi sources and write the closure report to `dev/gap_reports/gemmi_bio_remaining_full_port_closure.md` with every unsupported source field and every remaining non-green marker enumerated.
Step 72 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 73 [x]: Run `cargo fmt --all`.
Step 74 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 75 [x]: Run `cargo check -p cosmolkit-core --features op-contracts-strict`.
Step 76 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 77 [x]: Run `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
Step 78 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 79 [x]: Run `cargo test --workspace --release --features cosmolkit-core/op-contracts-strict`.
Step 80 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 81 [x]: Run `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html`.
Step 82 [x]: Read `dev/policy_invariants.md` and `dev/source_reproduction_protocol.md` to reload and follow the required execution standard, source reproduction rules, artifact requirements, no-git rule, and completion criteria for the next task.
Step 83 [x]: Run `.venv/bin/basedpyright python/tests python/examples`.
