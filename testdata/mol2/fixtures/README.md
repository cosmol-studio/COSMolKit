# MOL2 Fixtures

`rdkit/` contains byte-for-byte copies from RDKit `2026.03.1`, revision
`351f8f378f8ad6bbd517980c38896e66bf907af8c`. They are the complete MOL2
fixture set selected by the COSMolKit MOL2 parity suite. RDKit-origin files
retain RDKit's BSD-3-Clause license context.

`3rj7_ligand.mol2` is a pre-migration COSMolKit regression fixture whose
embedded name refers to PDB structure `3RJ7`. The file records its creating
user and 2021 export timestamps, but the exact download source and export
command were not preserved. It must not be described as a byte-for-byte RCSB
or RDKit copy until that history is recovered.

`source_manifest.jsonl` records this distinction and the byte length and
SHA-256 of every fixture.
