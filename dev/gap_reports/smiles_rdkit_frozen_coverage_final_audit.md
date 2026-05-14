# SMILES RDKit Frozen Coverage Final Audit

Date: 2026-05-14

## Verdict

Frozen-scope coverage is **not** exactly `100%`.

This final audit therefore does **not** certify checklist completion beyond the
frozen-scope reporting/documentation steps. The current tree still contains
remaining frozen-baseline copied-source gaps, and the final strict test run does
not pass cleanly.

## Basis For Verdict

### 1. Frozen-scope gap audit remains open

The frozen-baseline rescan documented in
`dev/gap_reports/smiles_rdkit_remaining_source_scan.md` found unresolved marker
classes inside the frozen scope, including:

- `notation/smiles.rs`: 1 `RDKit❌❌`, 2 `RDKit❗❗`, 14 `RDKit✔️❌`, 2037 `RDKit❗✔️`
- `notation/smiles_write.rs`: 2 `RDKit❗❗`, 582 `RDKit❗✔️`
- `notation/canon_rank.rs`: 119 `RDKit✔️❌`, 83 `RDKit❗✔️`
- `chemistry/rings.rs`: 7 `RDKit❌❌`, 78 `RDKit❗✔️`
- `chemistry/kekulize.rs`: 365 `RDKit❗✔️`
- `chemistry/valence.rs`: 13 `RDKit❗✔️`
- `operations/ops.rs`: 3 `RDKit❗❗`, 180 `RDKit✔️❌`, 33 `RDKit❗✔️`

Because the frozen scope still contains unresolved or known-gap copied-source
blocks, the result cannot be `100%`.

### 2. Final strict test run failed

`cargo test -p cosmolkit-core --features op-contracts-strict` completed with
failures in the final validation pass.

Observed failing tests:

- `io::sdf::tests::read_sdf_from_str_reads_simple_v3000_ctab`
- `notation::smiles::tests::assign_chiral_types_from_bond_dirs_promotes_implicit_h_on_three_coordinate_center`
- `notation::smiles::tests::handle_cx_part_and_name_clears_wiggly_single_bond_dirs_and_marks_unknown_stereo_like_rdkit`
- `notation::smiles::tests::handle_cx_part_and_name_parses_double_bond_stereo_like_rdkit`
- `notation::smiles::tests::mol_from_smiles_skip_cleanup_preserves_parser_temporary_props_like_rdkit_reader_flag`
- `notation::smiles_write::tests::cx_individual_atom_labels_writes_label_entries`
- `notation::smiles_write::tests::cx_individual_coords_writes_atom_order_coordinates`
- `notation::smiles_write::tests::isomeric_writer_outputs_default_non_tetrahedral_chiral_classes_like_rdkit`
- `notation::smiles_write::tests::isomeric_writer_outputs_non_tetrahedral_chiral_classes_like_rdkit`

The final strict test result was:

- 730 passed
- 9 failed

That alone prevents a `100% coverage` certification under the checklist's
closure standard.

## Files Updated In This Final Audit Segment

- `dev/gap_reports/smiles_rdkit_remaining_source_scan.md`
- `dev/porting_inventory.md`
- `crates/cosmolkit-core/src/support.rs`

## Final Statement

The frozen baseline has been audited and the repository documentation has been
updated to match the current implementation state, but the frozen-scope result
is **below** `100%` and must continue to be treated as incomplete.
