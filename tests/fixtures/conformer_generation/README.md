Copied RDKit conformer-generation fixtures live under `rdkit/test_data/`.

- Upstream library: RDKit `2026.03.1`
- Upstream source tree: `third_party/rdkit`
- Pinned upstream revision: `Release_2026_03_1` (`351f8f378f8ad6bbd517980c38896e66bf907af8`)
- Provenance manifest: `rdkit_inventory.jsonl`

The vendored files are byte-for-byte copies of
`third_party/rdkit/Code/GraphMol/DistGeomHelpers/test_data/*` so default Rust
tests and golden generation do not require the RDKit submodule to be checked
out in CI.
