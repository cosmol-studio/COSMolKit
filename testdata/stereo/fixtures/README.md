# Stereo Fixtures

This directory stores committed stereo input fixtures and provenance, not
generated expectations. Preparation, preflight, and comparison are managed by
[`parity-tests_fixed`](../../../parity-tests_fixed/README.md).

Recorded sources, byte lengths, and SHA-256 checksums are in
[`source_manifest.jsonl`](source_manifest.jsonl).

## Ligand reader regressions

[`1aid_ligand.mol2`](1aid_ligand.mol2) and
[`1aid_ligand.sdf`](1aid_ligand.sdf) preserve ligand stereochemistry across
MOL2 and SDF readers. They predate the migration and refer to PDB `1AID`;
their headers identify `X-TOOL` (2018) and `I-interpret`, respectively. The
download source and conversion commands are unknown: these are not verified
byte-for-byte upstream copies.

## AssignAtomChiralTagsFromStructure cases

[`assign_atom_chiral_tags_from_structure_cases.json`](assign_atom_chiral_tags_from_structure_cases.json)
contains 47 primary and 30 octahedral cases authored by COSMolKit against
RDKit `2026.03.1`, not copied upstream data. Source selection and observable
state are documented in the
[source audit](../../../dev/gap_reports/rdkit_assign_atom_chiral_tags_from_structure_source_audit.md).
Numeric strings specify exact binary64 construction, not tolerances;
octahedral cases cover every nested switch branch with both volume signs.
All cases, including exceptions and non-finite coordinates, must be retained.
Execution and generated identity checks belong to the
[special-regression lane](../../../parity-tests_fixed/README.md#special-regressions),
not a SMILES profile.

## Modern CIPLabeler cases

[`ciplabeler_focused.json`](ciplabeler_focused.json) is a COSMolKit-authored
matrix for RDKit `2026.03.1` modern `CIPLabeler::assignCIPLabels`, using pinned
upstream molecular inputs where available. It covers repeated calls, selection
masks, recursion profiles, unsupported dispatch, and all ten descriptors:
`R`, `S`, `r`, `s`, `E`, `Z`, `M`, `P`, `m`, `p`. Generated expectations
record complete observable properties, stereo state, and success/errors after
each call.

## Python stereoisomer enumeration cases

[`rdkit_python_stereoisomer_cases.json`](rdkit_python_stereoisomer_cases.json)
defines inputs/options for the pinned Python `FindPotentialStereo` → flippers
→ `EnumerateStereoisomers` boundary, with per-case source locators and branch
rationale. Expectations are generated separately; no mismatch is preaccepted.

[`rdkit/two_centers_or.mol`](rdkit/two_centers_or.mol) and
[`rdkit/simple_either.mol`](rdkit/simple_either.mol) are byte-for-byte pinned
RDKit FileParsers fixtures used by `UnitTestMol3D.py`. Tests use these committed
copies without a runtime dependency on `third_party/rdkit`.
