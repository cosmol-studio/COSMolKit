# Descriptor Test Data

`fixtures/rdkit/high_feasibility_descriptor_focused.smi` is the committed
source-profile input for the high-feasibility descriptor parity suite. Its
comments identify the covered behavior and the pinned RDKit regression source.
The explicit `<EMPTY>` row denotes the empty SMILES string; blank lines remain
formatting and comments remain provenance only.

This directory records fixture provenance. Current corpus preparation
and comparison use the [standard test runner](../../parity-tests/README.md),
not an owner-local generator. Generated references remain outside Git and
record the reference version, inputs, parameters, schema, counts and checksums.

The upstream anchors are RDKit 2026.03.1 revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`:

- `Code/GraphMol/Descriptors/test.cpp`
- `Code/GraphMol/Descriptors/Wrap/testMolDescriptors.py`
- `rdkit/Chem/UnitTestGraphDescriptors_2.py`
