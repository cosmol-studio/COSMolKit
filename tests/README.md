# Test Layout

- `fixtures/chem/` chemistry test molecules and expected outputs
- `fixtures/bio/` structure biology test inputs and expected outputs
- `fixtures/mol2/rdkit/` copied RDKit MOL2 parser fixtures used by MOL2 read parity; this avoids requiring CI to see `third_party/rdkit` as a checked-out submodule when regenerating goldens
- `fixtures/conformer_generation/rdkit/test_data/` copied RDKit conformer-generation fixtures used by conformer parity and golden regeneration; provenance is recorded in `fixtures/conformer_generation/rdkit_inventory.jsonl`
- `smiles.smi` shared SMILES corpus for graph-feature and molblock parity tests
- `golden/graph_features.jsonl` RDKit baseline for atom, bond, valence, stereo, and CIP graph-feature parity
- `golden/molblock_v2000_minimal.jsonl` RDKit baseline for minimal V2000 mol block body parity
- `golden/molblock_v2000_kekulized.jsonl` RDKit baseline for kekulized bond-block parity (ignores coordinates)
- `golden/tetrahedral_stereo_geometry.jsonl` RDKit ETKDG geometry baseline for tetrahedral stereo volume checks
- `golden/smiles_writer.jsonl` RDKit baseline for `MolToSmiles()` parity across `isomericSmiles`, `kekuleSmiles`, and `canonical` branches
- `golden/isomeric_smiles.jsonl` RDKit baseline for focused isomeric SMILES parity cases
- `golden/dg_bounds_matrix.jsonl` RDKit baseline for distance-geometry bounds matrix parity
- `golden/forcefield_params.jsonl` RDKit baseline for UFF/MMFF parameter coverage, initial energy/gradient, and single-/multi-conformer force-field optimization parity
- `golden/morgan_fingerprint.jsonl` RDKit baseline for Morgan fingerprint bit-vector and adjacent-row Tanimoto parity across radius, bit-count, chirality, bond-type, count-simulation, custom-invariant, generator, and AdditionalOutput branches
- `golden/svg_drawer.jsonl` RDKit baseline for `MolDraw2DSVG` output parity
- `golden/prepared_draw_molecule.jsonl` RDKit baseline for `PrepareMolForDrawing(kekulize=True, addChiralHs=True, wedgeBonds=True, forceCoords=True)` prepared atom coordinates and bond directions
- `golden/sdf_write.jsonl` RDKit baseline for MolBlock/SDF write parity across 2D/3D, V2000/V3000, stereo, and kekulize branches
- `golden/sdf_read.jsonl` RDKit baseline for SDF read parity across 2D/3D, V2000/V3000, stereo-marker, and coordinate-inferred branches
- `golden/molfile_read.jsonl` RDKit baseline for direct `.mol`/molblock read parity across the same CTAB branch matrix without SDF record separators or data fields
- `golden/mol2_read.jsonl` RDKit baseline for MOL2 read parity over copied RDKit MOL2 fixtures, covering parser parameters, topology, atom/bond fields, 3D coordinates, chirality, and SMILES output
- `golden/xyz_read.jsonl` RDKit baseline for XYZ block read parity: atom identities, one 3D conformer, coordinates, and no inferred bonds
- `corpus/topology/core.csv` target contract corpus for topology-changing operation invariant tests
- `corpus/topology/cow_small.csv` small COW-only topology corpus; do not run broad parity matrices just to prove value isolation
- `known_failures/topology_invariants.jsonl` exact xfail records for topology invariant failures; records must match operation, case, invariant, and error kind

## Standard Workflow (RDKit Parity)

RDKit `2026.03.1` is the current oracle for generated golden files. The source reference is `third_party/rdkit` pinned to `Release_2026_03_1` (`351f8f378f8ad6bbd517980c38896e66bf907af8`). Keep the Python environment project-level so the same `.venv` can later host COSMolKit Python bindings for direct comparison.

RDKit parity tests are strict source-level reproduction tests against `third_party/rdkit`. Do not make parity tests pass by loosening assertions, skipping mismatching fields, adding vague fallbacks, simplifying test conditions, row-specific patches, or heuristic guesses. When a mismatch appears, locate the corresponding RDKit source path and port that behavior directly; if the path is not implemented yet, keep the failure explicit and narrowly described.

1. Create project-level Python env (one shared env for testing + future bindings):
   - `uv sync --group dev`
   - this creates/updates `.venv/` at repository root
2. Regenerate all RDKit golden files through the single entrypoint:
   - `.venv/bin/python tests/scripts/gen_all_rdkit_goldens.py --python .venv/bin/python --clean --jobs 4`
   - omit `--jobs` to use the script default (`min(4, cpu_count, generator_count)`), or pass a larger value for a local high-core machine
3. Do not regenerate committed golden files through individual generator scripts.
   The per-surface scripts are implementation details of the unified entrypoint;
   using one-off outputs as committed baselines risks drift between parity
   surfaces.
4. (Optional) Install local COSMolKit Python build into the same env for direct comparison:
   - `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
5. Run core tests. Use a focused debug-profile filter while iterating, and
   release mode with the same strict feature set for full local parity runs:
   - `cargo test -p cosmolkit-core --features op-contracts-strict <test-filter>`
   - `cargo test -p cosmolkit-core --release --features op-contracts-strict`
6. Run the exhaustive SMILES writer branch matrix when checking full writer parity:
   - `cargo test -p cosmolkit-core --release --features op-contracts-strict --test rdkit_smiles_writer_parity smiles_writer_matches_rdkit_golden_across_param_branches -- --ignored`

Notes:
- `cargo test` never generates RDKit golden files. If a golden file is missing, regenerate all RDKit goldens with the single entrypoint before running Rust tests.
- CI regenerates all RDKit Python golden files through `tests/scripts/gen_all_rdkit_goldens.py` immediately after `uv sync --group dev` and before running Rust tests.
- The unified golden generator includes an RDKit version assertion for ETKDG geometry goldens so test conditions do not drift silently.

## Gemmi Macromolecular Parsing

Gemmi `v0.7.5` is the source-level reproduction target for PDB/mmCIF macromolecular parsing. The source reference is `third_party/gemmi` pinned to `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`. Gemmi comparison fixtures are explicit and separate from RDKit chemistry parity fixtures.

## Topology Operation Invariants

The topology invariant corpus and xfail records live under `tests/corpus/topology/` and `tests/known_failures/topology_invariants.jsonl`. They describe the target contract for topology-changing operations and are used by operation-system development work.

The target `_with_report()` operation design returns the transformed molecule plus atom/bond mappings and policy metadata. The invariant runner is intended to verify graph structure, conformer alignment, cache freshness, stereo reference validity, property preservation/remapping, copy-on-write isolation, optional I/O roundtrip, and optional RDKit parity.

COW is a sub-invariant, not a reason to run the whole RDKit parity corpus. Use `corpus/topology/cow_small.csv` for value-isolation coverage and `corpus/topology/core.csv` for the full topology operation contract.

Known topology failures are xfail records, not skips. The case still runs, the failure shape must match the record exactly, and an unexpected pass should fail so the record is removed.

`crates/cosmolkit-core/tests/rdkit_graph_feature_parity.rs` contains:
- `graph_feature_golden_has_one_record_per_smiles`
- `graph_feature_golden_preserves_rdkit_cip_fields_for_chiral_atoms`
- `graph_features_match_rdkit_golden_for_direct_and_explicit_hydrogen_molecules`

The graph feature test compares both direct molecules and explicit-hydrogen molecules. It covers atom atomic number, chirality, CIP code/rank, degree, formal charge, total hydrogens, radical electrons, hybridization, aromaticity, ring membership, and bond type/stereo/conjugation.

`crates/cosmolkit-core/src/io/molblock.rs` contains:
- `molblock_v2000_body_matches_rdkit_coordinates_and_topology`
- `molblock_kekulized_topology_matches_rdkit_golden`
- focused diagnostics for source-aligned depiction work

`crates/cosmolkit-core/tests/tetrahedral_stereo_geometry.rs` contains:
- `tetrahedral_stereo_ordered_ligands_match_rdkit_etkdg_positive_volume`
- missing-golden checks for the ETKDG geometry golden
- oriented-volume validation for `Molecule::tetrahedral_stereo()` (spec: `dev/tetrahedral_stereo.md`)

`crates/cosmolkit-core/tests/rdkit_dg_bounds_parity.rs` contains:
- `dg_bounds_golden_has_one_record_per_smiles`
- `dg_bounds_matrix_matches_rdkit_golden`
- strict RDKit parity coverage for distance-geometry bounds generation

`crates/cosmolkit-core/tests/rdkit_conformer_generation_library_parity.rs` contains:
- `conformer_generation_library_golden_has_one_record_per_smiles`
- `conformer_generation_library_matches_rdkit_golden_under_fixed_iteration_budget`
- corpus-wide RDKit golden coverage for the explicit path `SMILES -> AddHs -> ETKDGv3(seed=61453,numThreads=1,timeout=10) -> single conformer`

`crates/cosmolkit-core/tests/rdkit_forcefield_params_parity.rs` contains:
- `forcefield_params_golden_has_one_record_per_smiles_library_entry`
- `uff_has_all_molecule_params_matches_rdkit_golden`
- `mmff_has_all_molecule_params_matches_rdkit_golden`
- `mmff_atom_types_match_rdkit_golden`
- `forcefield_initial_energy_matches_rdkit_golden_for_first_embedded_row`
- `forcefield_initial_gradient_matches_rdkit_golden_for_first_embedded_row`
- `uff_single_conformer_final_coordinates_match_rdkit_golden_for_first_embedded_row`
- `mmff_single_conformer_final_coordinates_match_rdkit_golden_for_first_embedded_row`
- `uff_multi_conformer_final_coordinates_match_rdkit_golden_for_first_embedded_row`
- `mmff_multi_conformer_final_coordinates_match_rdkit_golden_for_first_embedded_row`
- strict RDKit parity coverage for UFF/MMFF parameter coverage, initial energy/gradient, and single-/multi-conformer optimization on the shared corpus

`crates/cosmolkit-core/src/properties/fingerprint.rs` contains focused Morgan
fingerprint unit coverage for deterministic bit generation, chirality,
count-simulation, custom atom/bond invariants, feature invariants, `fromAtoms`,
`ignoreAtoms`, Tanimoto, AdditionalOutput, and parameter validation. The
generated RDKit baseline `tests/golden/morgan_fingerprint.jsonl` is kept for
parity data, but the current core crate does not expose a dedicated
`rdkit_morgan_fingerprint_parity.rs` integration test. Python binding coverage
lives in `python/tests/test_morgan_fingerprint.py` and covers default Morgan
fingerprints, chirality, count simulation, feature/connectivity generators,
custom invariants, `fromAtoms`, batch fingerprint lists, Tanimoto similarity,
and AdditionalOutput fields.

`crates/cosmolkit-core/tests/rdkit_smiles_writer_parity.rs` contains:
- `smiles_writer_golden_has_one_record_per_smiles`
- `smiles_writer_matches_rdkit_golden_across_param_branches`
- strict RDKit parity coverage for `MolToSmiles()` across `isomericSmiles`, `kekuleSmiles`, and `canonical` branches

`crates/cosmolkit-core/tests/rdkit_svg_draw_parity.rs` contains:
- `svg_draw_golden_has_one_record_per_smiles`
- `svg_drawer_matches_rdkit_golden_except_tool_identifiers`
- `svg_drawer_matches_rdkit_golden_except_tool_identifiers_in_parallel_batch`
- strict RDKit parity coverage for the Python-exposed `MolDraw2DSVG(..., noFreetype=True)` final SVG string, excluding expected tool identifier metadata differences

`crates/cosmolkit-core/tests/rdkit_sdf_read_parity.rs` contains:
- `sdf_read_golden_covers_expected_case_matrix`
- `sdf_read_representative_apis_and_parameters_match_rdkit`
- `sdf_read_all_apis_and_parameters_match_rdkit` (`#[ignore]` exhaustive matrix)

`crates/cosmolkit-core/tests/rdkit_molfile_read_parity.rs` contains:
- `molfile_read_topology_and_atom_fields_match_rdkit`
- `molfile_read_coordinates_match_rdkit_for_2d_and_3d_records`
- `molfile_read_chirality_matches_rdkit_for_markers_and_coordinates`
- `molfile_read_to_smiles_matches_rdkit_canonical_and_noncanonical`

`crates/cosmolkit-core/tests/rdkit_mol2_read_parity.rs` contains:
- `mol2_read_golden_covers_expected_case_matrix`
- `mol2_read_topology_atom_and_bond_fields_match_rdkit`
- `mol2_read_coordinates_and_chirality_match_rdkit`
- `mol2_read_to_smiles_matches_rdkit_canonical_and_noncanonical`

PDB/mmCIF element-table and metal handling coverage includes:
- `rdkit_atomic_number_from_symbol_covers_canonical_and_legacy_entries`
- `bio_element_symbol_lookup_covers_rdkit_periodic_table`
- `pdb_element_field_accepts_full_rdkit_gemmi_element_table`
- `mmcif_atom_site_accepts_full_rdkit_gemmi_element_table`
- `rdkit_pdb_molecule_conversion_accepts_full_pdb_element_table`
- `pdb_connection_helpers_use_gemmi_full_metal_table`
- Python coverage for `MoleculeEdit.add_atom("Hg")`,
  `Molecule.from_pdb_block()` Hg/Cd records, and `Protein.from_pdb_str()` with
  supported metal HETATM records.

Current status:
- `cosmolkit-core` graph-feature parity is currently passing on the shared corpus (direct + explicit-H comparisons).
- tetrahedral stereo ordered-ligand geometry validation is currently passing against RDKit ETKDGv3 (`seed=42`) on all chiral corpus entries.
- DG bounds matrix parity is currently passing on the shared corpus.
- force-field parameter parity is currently passing on the shared corpus across UFF/MMFF parameter coverage, initial energy/gradient, and single-/multi-conformer optimization branches.
- Morgan fingerprint unit and Python binding coverage is currently passing for
  the supported public branches; a dedicated current-core RDKit golden parity
  integration test is not present.
- SMILES writer parity is currently passing on the shared corpus across `isomericSmiles`, `kekuleSmiles`, and `canonical` branches.
- strict V2000 molblock coordinate/topology parity is currently passing on the shared corpus.
- MOL2 read parity is currently passing against copied RDKit MOL2 parser fixtures for the generated MOL2 golden.
- SVG drawer parity is currently passing against RDKit final SVG goldens for the shared corpus.
- Temporary stress check result: random sampling 1000 SMILES from `core_comp_lib.csv` with regenerated RDKit goldens still exposes unresolved molblock parity gaps (details logged under `tmp/rust_test_core_comp_lib_sample1000_with_regen_errors.txt`).
