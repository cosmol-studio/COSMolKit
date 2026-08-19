# MolBlock Stress-Regressions Fixtures

`molblock_stereo_direction_components.jsonl` contains 14 RDKit-generated
V3000 MolBlocks: deterministic 2D and 3D forms for each of the 7 unique
structures whose parsed directional-bond state remained allocator-dependent
after source-level MolBlock fixes.

The source structures came from the ChEMBL 37 structure table used by the
2026-08-17 stress run. The complete 2,897,819-record source export had SHA-256
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.
Selection retained every unique source row whose fresh replay differed only in
absolute `ENDUPRIGHT`/`ENDDOWNRIGHT` representation. Fixture MolBlocks and
expected parsed state were produced with pinned RDKit `2026.03.1`, 2D
coordinate generation plus wedging, or ETKDGv3 with random seed 42 and one
thread, followed by V3000 writing with stereo and kekulization enabled and an
RDKit sanitize/reparse with explicit hydrogens retained.

RDKit stores tied directional bonds in an ordered container whose pointer
tie-break depends on allocator history. The regression therefore requires all
stable atom, bond, E/Z, stereo-control, aromaticity, canonical-SMILES, and
coordinate fields to match while permitting only one uniform up/down inversion
per connected stereogenic-double-bond constraint component. A local or
unconstrained inversion is not accepted.

The fixture is derived from ChEMBL data and retains that data's license
context; it is not copied from RDKit source code. `source_manifest.jsonl`
records its exact checksum and selection identity.
