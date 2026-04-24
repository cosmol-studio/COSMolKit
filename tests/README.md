# Test Layout

- `fixtures/chem/` chemistry test molecules and expected outputs
- `fixtures/bio/` structure biology test inputs and expected outputs
- `smiles.txt` SMILES corpus for parity tests
- `golden/atomic_nums.jsonl` RDKit baseline for `atom.GetAtomicNum()` parity
- `molblock_minimal_smiles.txt` corpus for minimal V2000 body parity
- `golden/molblock_v2000_minimal.jsonl` RDKit baseline for minimal V2000 mol block body parity
- `golden/molblock_v2000_kekulized.jsonl` RDKit baseline for kekulized bond-block parity (ignores coordinates)

## Standard Workflow (SMILES Atomic Number Parity)

1. Create project-level Python env (one shared env for testing + future bindings):
   - `uv sync --group dev`
   - this creates/updates `.venv/` at repository root
2. Regenerate RDKit golden:
   - `.venv/bin/python tests/scripts/gen_rdkit_atomic_nums.py --input tests/smiles.txt --output tests/golden/atomic_nums.jsonl`
   - `.venv/bin/python tests/scripts/gen_rdkit_v2000_minimal_golden.py --input tests/molblock_minimal_smiles.txt --output tests/golden/molblock_v2000_minimal.jsonl`
   - `.venv/bin/python tests/scripts/gen_rdkit_kekulize_molblock_golden.py --input tests/molblock_minimal_smiles.txt --output tests/golden/molblock_v2000_kekulized.jsonl`
3. (Optional) Install local COSMolKit Python build into the same env for direct comparison:
   - `.venv/bin/maturin develop --manifest-path bindings/python/native/Cargo.toml`
4. Run Rust tests:
   - `cargo test -p cosmolkit-chem-core`
5. Run one-off feature parity check for a single SMILES (ours vs RDKit):
   - `cargo run -p cosmolkit-chem-core --example feature_parity -- "C[C@H](F)Cl"`

Notes:
- `cargo test -p cosmolkit-chem-core` will auto-generate `tests/golden/atomic_nums.jsonl` if it is missing.
- `cargo test -p cosmolkit-io` will auto-generate `tests/golden/molblock_v2000_minimal.jsonl` if it is missing.
- `cargo test -p cosmolkit-io` will auto-generate `tests/golden/molblock_v2000_kekulized.jsonl` if it is missing.
- Python lookup order for auto-generation: `COSMOLKIT_PYTHON` -> `.venv/bin/python` -> `python3`.

`crates/chem-core/tests/rdkit_atomic_num_parity.rs` contains:
- `golden_has_entry_for_each_smiles_input` (always enabled, CI-safe)
- `atomic_numbers_match_rdkit_golden` (currently `#[ignore]`; enable after SMILES parser implementation)
