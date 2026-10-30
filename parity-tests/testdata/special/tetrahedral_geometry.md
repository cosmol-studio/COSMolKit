# Tetrahedral geometry inputs

The 30 SMILES come from the historical 31-row `smiles_small`
`tetrahedral_stereo_geometry.jsonl`. At the user's request, `historical_23`
is excluded: RDKit ETKDGv3 with seed 42 returned -1, leaving no reference
geometry to compare. Remaining case IDs are unchanged. The failed run is
recorded in `parity-tests/reports/run-iWfejj`; no other inputs are filtered.
Source snapshot: RDKit 2026.03.1, ETKDGv3, random seed 42;
SHA-256 `6df02f48431d190afc410640af85129e6ed318e9ae41a52cd2cee3d53da952a0`.
Only SMILES and preparation options are retained here, not historical outputs.
The regular prepare command generates references using the runner's current pin.

The predicate and scalar/batch coverage restore commit `24cfdb1`'s
`tetrahedral_stereo_ordered_ligands_match_rdkit_etkdg_positive_volume`.
See the runner README for the shared prepare/compare workflow.
