# Test Layout

- `fixtures/chem/` chemistry test molecules and expected outputs
- `fixtures/bio/` structure biology test inputs and expected outputs
- `smiles.smi` shared SMILES corpus for graph-feature and molblock parity tests
- `golden/graph_features.jsonl` RDKit baseline for atom, bond, valence, stereo, and CIP graph-feature parity
- `golden/molblock_v2000_minimal.jsonl` RDKit baseline for minimal V2000 mol block body parity
- `golden/molblock_v2000_kekulized.jsonl` RDKit baseline for kekulized bond-block parity (ignores coordinates)
- `golden/tetrahedral_stereo_geometry.jsonl` auto-generated RDKit ETKDG geometry baseline for tetrahedral stereo volume checks
- `golden/smiles_writer.jsonl` RDKit baseline for `MolToSmiles()` parity across `isomericSmiles`, `kekuleSmiles`, and `canonical` branches
- `golden/dg_bounds_matrix.jsonl` RDKit baseline for distance-geometry bounds matrix parity

## Standard Workflow (RDKit Parity)

RDKit `2026.03.1` is the current oracle for generated golden files. The source reference is `third_party/rdkit` pinned to `Release_2026_03_1` (`351f8f378f8ad6bbd517980c38896e66bf907af8`). Keep the Python environment project-level so the same `.venv` can later host COSMolKit Python bindings for direct comparison.

RDKit parity tests are strict source-level reproduction tests against `third_party/rdkit`. Do not make parity tests pass by loosening assertions, skipping mismatching fields, adding vague fallbacks, simplifying test conditions, row-specific patches, or heuristic guesses. When a mismatch appears, locate the corresponding RDKit source path and port that behavior directly; if the path is not implemented yet, keep the failure explicit and narrowly described.

1. Create project-level Python env (one shared env for testing + future bindings):
   - `uv sync --group dev`
   - this creates/updates `.venv/` at repository root
2. Regenerate RDKit golden:
   - `.venv/bin/python tests/scripts/gen_rdkit_graph_features.py --input tests/smiles.smi --output tests/golden/graph_features.jsonl`
   - `.venv/bin/python tests/scripts/gen_rdkit_v2000_minimal_golden.py --input tests/smiles.smi --output tests/golden/molblock_v2000_minimal.jsonl`
   - `.venv/bin/python tests/scripts/gen_rdkit_kekulize_molblock_golden.py --input tests/smiles.smi --output tests/golden/molblock_v2000_kekulized.jsonl`
   - `.venv/bin/python tests/scripts/gen_rdkit_tetrahedral_stereo_geometry.py --input tests/smiles.smi --output tests/golden/tetrahedral_stereo_geometry.jsonl`
   - `.venv/bin/python tests/scripts/gen_rdkit_smiles_writer_golden.py --input tests/smiles.smi --output tests/golden/smiles_writer.jsonl`
   - `.venv/bin/python tests/scripts/gen_rdkit_dg_bounds_golden.py --input tests/smiles.smi --output tests/golden/dg_bounds_matrix.jsonl`
3. (Optional) Install local COSMolKit Python build into the same env for direct comparison:
   - `.venv/bin/maturin develop --manifest-path python/Cargo.toml`
4. Run Rust tests:
   - `cargo test -p cosmolkit-core`
   - `cargo test -p cosmolkit`
5. Run one-off feature parity check for a single SMILES (ours vs RDKit):
   - `cargo run -p cosmolkit-core --example feature_parity -- "C[C@H](F)Cl"`
6. Generate one-off minimal SDF output:
   - `cargo run -p cosmolkit-core --example minimal_sdf -- "CCO"`

Notes:
- `cargo test -p cosmolkit-core` will auto-generate `tests/golden/graph_features.jsonl` if it is missing.
- `cargo test -p cosmolkit-core` will auto-generate `tests/golden/molblock_v2000_minimal.jsonl` if it is missing.
- `cargo test -p cosmolkit-core` will auto-generate `tests/golden/molblock_v2000_kekulized.jsonl` if it is missing.
- `cargo test -p cosmolkit-core --test tetrahedral_stereo_geometry` will auto-generate `tests/golden/tetrahedral_stereo_geometry.jsonl` if it is missing.
- `cargo test -p cosmolkit-core --test rdkit_smiles_writer_parity` will auto-generate `tests/golden/smiles_writer.jsonl` if it is missing.
- `cargo test -p cosmolkit-core --test rdkit_dg_bounds_parity` will auto-generate `tests/golden/dg_bounds_matrix.jsonl` if it is missing.
- Python lookup order for auto-generation: `COSMOLKIT_PYTHON` -> `.venv/bin/python` -> `python3`.
- `tests/scripts/gen_rdkit_tetrahedral_stereo_geometry.py` asserts `rdkit == 2026.3.1` before generating ETKDG geometry golden so test conditions do not drift silently.

`crates/cosmolkit-core/tests/rdkit_graph_feature_parity.rs` contains:
- `graph_feature_golden_has_one_record_per_smiles`
- `graph_feature_golden_records_cip_for_chiral_atoms`
- `graph_features_match_rdkit_golden_for_direct_and_explicit_hydrogen_molecules`

The graph feature test compares both direct molecules and explicit-hydrogen molecules. It covers atom atomic number, chirality, CIP code/rank, degree, formal charge, total hydrogens, radical electrons, hybridization, aromaticity, ring membership, and bond type/stereo/conjugation.

`crates/cosmolkit-core/src/io/molblock.rs` contains:
- `molblock_v2000_body_matches_rdkit_coordinates_and_topology`
- `molblock_kekulized_topology_matches_rdkit_golden`
- focused diagnostics for source-aligned depiction work

`crates/cosmolkit-core/tests/tetrahedral_stereo_geometry.rs` contains:
- `tetrahedral_stereo_ordered_ligands_match_rdkit_etkdg_positive_volume`
- auto-generation hook for the ETKDG geometry golden
- oriented-volume validation for `Molecule::tetrahedral_stereo()` (spec: `tetrahedral_stereo_representation.md`)

`crates/cosmolkit-core/tests/rdkit_dg_bounds_parity.rs` contains:
- `dg_bounds_golden_has_one_record_per_smiles`
- `dg_bounds_matrix_matches_rdkit_golden`
- strict RDKit parity coverage for distance-geometry bounds generation

`crates/cosmolkit-core/tests/rdkit_smiles_writer_parity.rs` contains:
- `smiles_writer_golden_has_one_record_per_smiles`
- `smiles_writer_matches_rdkit_golden_across_param_branches`
- strict RDKit parity coverage for `MolToSmiles()` across `isomericSmiles`, `kekuleSmiles`, and `canonical` branches

Current status:
- `cosmolkit-core` graph-feature parity is currently passing on the shared corpus (direct + explicit-H comparisons).
- tetrahedral stereo ordered-ligand geometry validation is currently passing against RDKit ETKDGv3 (`seed=42`) on all chiral corpus entries.
- DG bounds matrix parity is currently passing on the shared corpus.
- SMILES writer parity is currently passing on the shared corpus across `isomericSmiles`, `kekuleSmiles`, and `canonical` branches.
- strict V2000 molblock coordinate/topology parity is currently passing on the shared corpus.
- `cargo test -p cosmolkit-core` is currently green: `26 passed`, `0 failed`, `0 ignored`.
- `cargo test -p cosmolkit` is currently green.
- Temporary stress check result: random sampling 1000 SMILES from `core_comp_lib.csv` with regenerated RDKit goldens still exposes unresolved molblock parity gaps (details logged under `tmp/rust_test_core_comp_lib_sample1000_with_regen_errors.txt`).
