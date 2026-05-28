# Gemmi PDB/mmCIF/mmJSON Final Audit

Audit target: `crates/cosmolkit-core/src/io/bio.rs`

Audit sources:

- `third_party/gemmi/src/pdb.cpp`
- `third_party/gemmi/include/gemmi/pdb.hpp`
- `third_party/gemmi/src/mmcif.cpp`
- `third_party/gemmi/include/gemmi/mmcif.hpp`
- `third_party/gemmi/include/gemmi/mmread.hpp`
- `third_party/gemmi/include/gemmi/mmcif_impl.hpp`
- `third_party/gemmi/src/polyheur.cpp`

## Summary

The Gemmi-aligned structural reader closure requested by
`dev/pdb_mmcif_gemmi_full_port_checklist.md` is now exposed through:

- `BioStructure::from_structure_str(...)`
- `BioStructure::from_str_with_format(...)`
- `BioStructure::from_pdb_str(...)`
- `BioStructure::from_pdb_str_with_params(...)`
- `read_mmcif_atom_site_subset_from_str(...)`

Fixture-backed integration coverage now exercises full-path PDB, mmCIF, and
chem-comp dispatch reads through the public `io::bio` surface. Targeted unit
coverage now also exercises the mmJSON dispatch path and the mmJSON-to-CIF
document conversion surface used by Gemmi `mmread.hpp` / `json.cpp`.

## Remaining non-`Gemmi✔️✔️` markers

`crates/cosmolkit-core/src/io/bio.rs` currently contains:

- `3` `Gemmi❌❌` markers
- `10` `Gemmi❗✔️` markers
- `0` `Gemmi❗❗` markers
- `0` `Gemmi✔️❌` markers

These are the remaining visible gaps after checklist Step 255.

## Remaining `Gemmi❌❌`

### PDB stream parser control flow

Location: `BioStructure::from_pdb_str_with_params(...)`

- duplicate implicit-model rejection for `ATOM/HETATM between models`

These markers are now isolated to the exact `st.find_model(num)` duplicate
implicit-model guard inside the redesigned flat-row builder path. The visible
behavior for duplicate implicit-model collisions is covered by targeted tests,
but the copied lines remain marked because the current Rust path still uses the
redesigned builder state rather than a source-identical object lookup at that
exact statement boundary.

## Remaining `Gemmi❗✔️`

### PDB finalization

Location: `PdbBioBuilder::finish(...)`

- `st.setup_cell_images();`

### PDB stream parser structural control flow

Location: `BioStructure::from_pdb_str_with_params(...)`

- implicit model creation path
- chain creation/reset path
- residue map lookup/creation path
- fallback element inference branch
- atom push path

These markers reflect that the Rust row-builder reproduces the modeled behavior
through redesigned flat-row structures rather than source-identical object
mutation.

## Conclusion

The public Gemmi-aligned PDB/mmCIF/mmJSON structural reader surface requested by the
checklist is implemented, documented, and covered by targeted unit tests plus
fixture-backed integration tests. The remaining non-green markers are visible
and bounded:

- mmJSON dispatch and mmJSON-backed structure reads are now implemented through
  a reproduced `read_mmjson_insitu(...)`-equivalent document conversion path.
- the only remaining `Gemmi❌❌` markers are the duplicated implicit-model guard
  lines inside the PDB stream parser.
- several PDB helper lines remain `Gemmi❗✔️` because the Rust redesign does not
  yet mirror Gemmi's exact object-level mutation flow at those copied lines.
- no remaining markers are hidden behind heuristic fallback behavior.
