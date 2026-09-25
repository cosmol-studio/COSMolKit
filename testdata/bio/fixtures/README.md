# Structural-Biology Fixtures

These three files are repository-authored synthetic integration fixtures.
They are deliberately small records assembled to exercise PDB, mmCIF, and
chemical-component fields through the public structure readers. They are not
copies from Gemmi or from real PDB entries: identifiers such as `1ABC` and
`9XYZ` are synthetic test labels.

The files were preserved byte-for-byte from the pre-migration COSMolKit test
suite and use the repository's license context. `source_manifest.jsonl`
records each byte length and SHA-256.

The separate [`gemmi_residues/`](gemmi_residues/README.md) directory contains
an upstream-derived fixed residue-table snapshot with its own provenance,
checksums and license. It is not one of the three synthetic structures above.
