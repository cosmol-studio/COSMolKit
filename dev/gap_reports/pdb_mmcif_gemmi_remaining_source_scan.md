# Gemmi PDB/mmCIF Remaining Source Scan

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

`io::bio` currently exposes two subset readers:

- `read_pdb_coordinate_subset_from_str(...)`
- `read_mmcif_atom_site_subset_from_str(...)`

They already reproduce part of Gemmi's PDB atom/ANISOU, SEQRES, HEADER/TITLE/KEYWDS/EXPDTA/AUTHOR/CRYST1, and mmCIF `_atom_site` / `_entity*` / `_struct_asym` behavior. They do not yet close the full Gemmi structural reader surface. The remaining work is still substantial and is concentrated in:

- exact PDB helper and remark parsing
- exact PDB record-branch closure and parser finalization
- mmCIF helper, metadata, connectivity, assembly, and chem-comp readers
- full-document dispatch and format detection helpers

## Current function-level status

### Already source-backed enough to keep, but still not full-reader closure

| Gemmi function / behavior | Current `io::bio` status | Notes |
|---|---|---|
| `populate_structure_from_pdb_stream(...)` ATOM/HETATM core field extraction | Partial | Basic atom rows, element inference, occupancy/B-factor, and residue creation exist; chain splitting, segment IDs, `after_ter`, and duplicate-model semantics remain open. |
| `populate_structure_from_pdb_stream(...)` ANISOU adjacency and tensor decode | Partial | Basic adjacency and duplicate checks are present; needs exact surrounding control flow closure with model handling. |
| `populate_structure_from_pdb_stream(...)` `SEQRES` | Partial | Entity sequence capture exists; later HELIX/SHEET/TER/finalization interactions are still missing. |
| `populate_structure_from_pdb_stream(...)` `HEADER` / `TITLE` / `KEYWDS` / `EXPDTA` / `AUTHOR` / `CRYST1` | Partial | Basic metadata fields and CRYST1 cell/H-M/Z_PDB parsing exist; matrix branches and exact full control flow remain open. |
| `read_atom_sites(...)` core column extraction | Partial | Basic atom-site rows, auth/label fallback for atom/comp/asym/seq, model grouping, occupancy/B-factor/formal charge are present. |
| `read_entity_and_sequence_info(...)` `_entity`, `_entity_poly`, `_entity_poly_seq`, `_struct_asym` | Partial | Polymer kind inference, microheterogeneity concatenation, and `_struct_asym` chain linking exist; `_struct_ref*` and no-`_struct_asym` fallback are still missing. |

### Remaining PDB helper functions

| Gemmi function | Remaining gap | Planned step |
|---|---|---|
| `gemmi::impl::set_cell_from_mmcif(...)` | Not ported into shared mmCIF top-level path. | 5 |
| `gemmi::impl::find_spacegroup_hm_value(...)` | Not ported into shared mmCIF top-level path. | 5 |
| `read_matrix(...)` | Not ported; exact row validation, short-line return, and matrix/vector layout are missing. | 11 |
| `add_software(...)` | Not ported. | 17 |
| `add_restraint_count_weight(...)` | Not ported. | 17 |
| `read_remark3_line(...)` | Not ported. | 17 |
| `read_remark_200_230_240(...)` | Not ported. | 17 |
| `complete_ssbond(...)` and helper atom completion behavior | Not ported. | 23 |
| `process_conn(...)` | Not ported. | 23 |
| `change_author_name_format_to_mmcif(...)` | Not ported. | 23 |
| `read_metadata_from_remarks(...)` | Not ported. | 29 |

### Remaining exact PDB parser control flow and record branches

| Gemmi function / branch group | Remaining gap | Planned step |
|---|---|---|
| `populate_structure_from_pdb_stream(...)` options and input guards | `max_line_length`, `split_chain_on_ter`, `skip_remarks`, `check_non_ascii`, CIF/mmJSON misdetection are not closed. | 35 |
| `populate_structure_from_pdb_stream(...)` model and atom-state control flow | Implicit model creation parity, duplicate-model rejection, `MODEL`/`ENDMDL`/`END`, and exact ANISOU state transitions are incomplete. | 41 |
| `populate_structure_from_pdb_stream(...)` chain and residue control flow | Chain splitting, segment IDs, residue-map reuse, and `after_ter` semantics are incomplete. | 47 |
| `populate_structure_from_pdb_stream(...)` `REMARK` and `CONECT` | Raw remark retention and deferred connection staging are absent. | 53 |
| `populate_structure_from_pdb_stream(...)` `HELIX` and `SHEET` | Not ported. | 59 |
| `populate_structure_from_pdb_stream(...)` `SSBOND`, `LINK`, `CISPEP` | Deferred record capture is absent. | 65 |
| `populate_structure_from_pdb_stream(...)` `TER`, `MODRES`, `HETNAM`, `DBREF*` | Not ported. | 71 |
| `populate_structure_from_pdb_stream(...)` `SCALEn`, `ORIGX`, `MTRIXn` | Not ported. | 77 |
| `assign_subchains(...)` parser-owned finalization path | Not ported. | 83 |
| `populate_structure_from_pdb_stream(...)` finalization tail, phase 1 | Empty-model insertion, TER cleanup, subchain propagation, and cell-image setup are absent. | 89 |
| `populate_structure_from_pdb_stream(...)` finalization tail, phase 2 | Deferred `process_conn`, author normalization, remark metadata import, and full CCD restoration are absent. | 95 |

### Remaining mmCIF helper and metadata functions

| Gemmi function | Remaining gap | Planned step |
|---|---|---|
| `copy_int(...)`, `copy_double(...)`, `copy_string(...)`, `get_smat33(...)`, `transform_tags(...)`, `get_transform_matrix(...)` | Not ported. | 101 |
| `make_seqid(...)`, `make_resid(...)`, `RowAccess`, `get_by_id(...)` | Not ported as exact helpers; local ad hoc fallback logic exists instead. | 107 |
| `get_anisotropic_u(...)` | Not ported. | 113 |
| `read_helices(...)` | Not ported. | 113 |
| `read_sheets(...)` | Not ported. | 113 |
| `set_part_of_address_from_label(...)` | Not ported. | 119 |
| `read_connectivity(...)` | Not ported. | 119 |
| `read_prot_cis(...)` | Not ported. | 119 |
| `read_struct_mod_residue(...)` | Not ported. | 119 |
| `parse_operation_expr(...)` | Not ported. | 125 |
| `read_assemblies(...)` | Not ported. | 125 |
| `fill_residue_entity_type(...)` | Not ported. | 125 |
| `read_sifts_unp(...)` | Not ported. | 125 |
| `find_diffrn(...)` | Not ported. | 125 |
| `read_entry_info(...)` | Not ported. | 131 |
| `read_audit_author(...)` | Not ported. | 131 |
| `read_refinement_info(...)` | Not ported. | 131 |
| `read_tls_info(...)` | Not ported. | 131 |
| `read_experimental_info(...)` | Not ported. | 131 |
| `read_reflns_info(...)` | Not ported. | 131 |
| `read_software_info(...)` | Not ported. | 131 |
| `read_ncs_info(...)` | Not ported. | 131 |

### Remaining mmCIF atom/entity reader gaps

| Gemmi function | Remaining gap | Planned step |
|---|---|---|
| `read_atom_sites(...)` | `label_entity_id`, anisotropic-U attachment, `group_PDB`, `calc_flag`, `pdbx_tls_group_id`, `ccp4_deuterium_fraction`, and inconsistent-sequence failure are missing. | 137 |
| `read_entity_and_sequence_info(...)` | `_struct_ref`, `_struct_ref_seq`, duplicate-dbref suppression, and `_struct_asym`-absence fallback are missing. | 143 |

### Remaining mmCIF top-level structure and chem-comp readers

| Gemmi function | Remaining gap | Planned step |
|---|---|---|
| `populate_structure_from_block(...)` | Top-level metadata reader ordering, fract transform, origx import, postprocessing, assembly/SIFTS import, and shortened CCD restoration are missing. | 149, 155 |
| `make_structure(cif::Document&&, cif::Document*)` | Multi-block `_atom_site` rejection and first-block-only coordinate semantics are missing. | 161 |
| `make_residue_from_chemcomp_block(...)` | Not ported. | 167 |
| `make_model_from_chemcomp_block(...)` | Not ported. | 167 |
| `make_structure_from_chemcomp_block(...)` | Not ported. | 173 |
| `check_chemcomp_block_number(...)` | Not ported. | 173 |
| `make_structure_from_chemcomp_doc(...)` | Not ported. | 173 |

### Remaining dispatch helpers from `mmread.hpp`

| Gemmi function | Remaining gap | Planned step |
|---|---|---|
| `coor_format_from_ext(...)` | Not ported. | 179 |
| `coor_format_from_content(...)` | Not ported as exact helper; current PDB reader does not expose the shared dispatch semantics. | 179 |
| `make_structure_from_doc(...)` | Not ported. | 185 |
| `read_structure_from_memory(...)` | Not ported. | 185 |

## Immediate implementation order confirmed by source scan

The next executable source-backed units are still the checklist front:

1. shared mmCIF cell / H-M extraction helpers
2. PDB matrix reader helper
3. PDB REMARK 3 / 200 / 230 / 240 helpers
4. PDB connection and REMARK 290 helpers
5. PDB remark metadata aggregation

That order is consistent with the current checklist and with the largest unresolved call sites already marked `Gemmi❌❌` in `io::bio`.
