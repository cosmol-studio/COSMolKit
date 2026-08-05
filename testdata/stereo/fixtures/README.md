# Stereo Fixtures

The `1aid_ligand.mol2` and `1aid_ligand.sdf` files are pre-migration
COSMolKit regression fixtures for preserving the same ligand stereochemistry
across MOL2 and SDF readers. Their names refer to PDB structure `1AID`; the
MOL2 file identifies `X-TOOL` and a 2018 creation date, and the SDF file
identifies `I-interpret`. The exact download source and conversion commands
were not preserved.

These files must not be presented as byte-for-byte upstream copies until that
history is recovered. `source_manifest.jsonl` records their current byte
lengths, SHA-256 values, and the explicit incomplete-provenance status.

## AssignAtomChiralTagsFromStructure cases

`assign_atom_chiral_tags_from_structure_cases.json` is a COSMolKit-authored,
committed input fixture for exact comparison with RDKit `2026.3.1`. It is not
copied from RDKit. The case selection comes from the active source closure in
`GraphMol/Chirality.cpp:45-54,860-900,3114-3500`,
`Geometry/point.h:147-234`, RDKit's focused chirality regressions, and
`dev/gap_reports/rdkit_assign_atom_chiral_tags_from_structure_source_audit.md`.

The fixture defines topology, bond insertion order, coordinates, conformer
ids/dimensionality, preexisting atom and molecule properties, public options,
and environment state. Numeric strings name exact binary64 construction rules
that the reference generator must decode; they are not tolerances. The
`octahedral_switch_cases` table covers every nested
`OctahedralPermFrom3D` switch branch with both volume signs. Expected outputs
are generated only through:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python --profile smiles_small --suite stereo
```

The generated rows must include all observable fields listed in the source
audit. No case may be filtered after RDKit execution, including exception and
non-finite-coordinate rows.
