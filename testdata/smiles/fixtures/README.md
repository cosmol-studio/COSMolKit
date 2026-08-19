# SMILES Stress-Regressions Fixtures

`rdkit_fragment_scope_stress_regressions.jsonl` contains the complete retained
set of canonical-SMILES rows that differed during the 2026-08-17 ChEMBL 37
stress run. It includes all 9 unique source structures and all 108 affected
writer profiles; no retained row was filtered after comparison.

The source structures came from the ChEMBL 37 structure table used by the
stress corpus. The complete 2,897,819-record source export had SHA-256
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.
Rows were selected only when pinned RDKit `2026.03.1` and COSMolKit `0.2.12`
produced different `MolToSmiles()` results in the exhaustive 8-boolean-option
and 3-root profile matrix. Expected strings are the retained RDKit results.

The selection was made by `tmp/parity-audit/audit_surfaces.py` as recorded in
the stress-run manifest. The original audit-script SHA-256 was
`39ca8e4bf2e251868f24548a938d5d4215d2cd906111aae780a5d5829513a788`.
The fixture is a focused regression input under the ChEMBL data license; it is
not copied from RDKit source code. `source_manifest.jsonl` records its checksum
and exact selection identity.
