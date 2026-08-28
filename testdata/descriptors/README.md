# Descriptor Test Data

`fixtures/rdkit/high_feasibility_descriptor_focused.smi` is the committed
source-profile input for the high-feasibility descriptor parity suite. Its
comments identify the covered behavior and the pinned RDKit regression source.
The explicit `<EMPTY>` row denotes the empty SMILES string; blank lines remain
formatting and comments remain provenance only.

Expected JSONL is generated through the repository-wide entrypoint:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python \
  --profile descriptors_focused \
  --suite descriptors
```

The generated `expected/rdkit/descriptors_focused/manifest.json` records the
pinned RDKit version, source revision, generator and input checksums, runtime,
schema, output checksum, and record count. Generated expected data remain
uncommitted under the repository-wide expected-data ignore rule.

The upstream anchors are RDKit 2026.03.1 revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`:

- `Code/GraphMol/Descriptors/test.cpp`
- `Code/GraphMol/Descriptors/Wrap/testMolDescriptors.py`
- `rdkit/Chem/UnitTestGraphDescriptors_2.py`
