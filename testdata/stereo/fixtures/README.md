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

## Modern CIPLabeler cases

`ciplabeler_focused.json` is a COSMolKit-authored case matrix for the pinned
RDKit `2026.03.1` modern `CIPLabeler::assignCIPLabels` boundary. Its molecular
inputs are taken from pinned upstream CIPLabeler regressions where available;
the selection, repeated-call sequences, selected atom/bond masks, recursion
profiles, and unsupported-dispatcher boundary are assembled by COSMolKit.

The generated oracle records the complete observable `_CIPComputed`,
`_CIPCode`, `_CIPNeighborOrder`, `_CIPRank`, chiral-tag, bond-stereo, stereo
atom, computed-property, success, and error state after every call. The focused
matrix covers all ten public descriptors emitted by the supported modern
assignment configurations: `R`, `S`, `r`, `s`, `E`, `Z`, `M`, `P`, `m`, and
`p`.

## Python stereoisomer enumeration cases

`rdkit_python_stereoisomer_cases.json` defines the focused input and option
matrix for the pinned Python `FindPotentialStereo` to flippers to
`EnumerateStereoisomers` boundary. Each case records its exact upstream source
locator, observable branch rationale, and applicable option profiles. Expected
records are generated separately; the input matrix contains no accepted
mismatch or expected-output shortcut.

`rdkit/two_centers_or.mol` and `rdkit/simple_either.mol` are byte-for-byte
copies of the pinned RDKit FileParsers fixtures used by `UnitTestMol3D.py` for
enhanced stereo and explicit either-double enumeration. Tests resolve the
committed copies and do not require `third_party/rdkit` at runtime.
