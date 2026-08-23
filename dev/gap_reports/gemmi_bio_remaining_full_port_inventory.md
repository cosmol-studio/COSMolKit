# Gemmi BioStructure Remaining Full-Port Inventory

## Baseline And Boundary

This audit compares the current `BioStructure` implementation with Gemmi
commit `5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`.

Primary sources:

- `third_party/gemmi/src/to_mmcif.cpp`
- `third_party/gemmi/include/gemmi/to_mmcif.hpp`
- `third_party/gemmi/include/gemmi/to_cif.hpp`
- `third_party/gemmi/include/gemmi/cifdoc.hpp`
- reader helpers already cited in `crates/cosmolkit-core/src/io/bio.rs`

The public ownership boundary remains unchanged: complete structural reading
and writing belong to `BioStructure`; `Protein` is a lossy amino-acid
projection and `Molecule` is a chemical graph. No structural writer belongs on
`Protein` or `Molecule`.

## Public Empty Surface

The only deliberately absent public structural-I/O surface is the mmCIF
writer family:

- `BioStructure::to_mmcif()`
- `BioStructure::to_mmcif_with_options(MmcifWriteOptions)`
- `BioStructure::write_mmcif(path)`
- `BioStructure::write_mmcif_with_options(path, MmcifWriteOptions)`

These methods cannot be added independently. Their complete source call path
is:

```text
MmcifOutputGroups
  -> make_mmcif_document
  -> update_mmcif_block
     -> add_cif_atoms
     -> write_assemblies
     -> write_ncs_oper
     -> write_struct_conn
     -> write_cispeps
  -> cif::write_cif_to_stream
     -> write_cif_block_to_stream
     -> write_out_item
     -> write_out_pair / write_out_loop / write_text_field
```

The current private `CifDocument`, `CifBlock`, `CifItem`, and `CifLoop` types
are reader-only containers. They lack Gemmi category replacement, pair/loop
interleaving, row insertion, quoting, and serialization behavior. A public
writer placed over them today would be a placeholder.

## `MmcifOutputGroups` Inventory

Gemmi exposes 34 output switches. All must exist with the source defaults:

| Group | Current state source | Port status |
|---|---|---|
| `atoms`, `group_pdb`, `auth_all` | atoms, residues, chains, models, coordinates | Missing writer |
| `block_name`, `entry` | `name`, `metadata.entry_id` | Missing writer |
| `database_status` | `metadata.received_initial_deposition_date` | Missing writer |
| `author` | `metadata.authors` | Missing writer |
| `cell`, `symmetry`, `scale` | `crystal` | Missing writer |
| `entity`, `entity_poly`, `entity_poly_seq`, `struct_asym` | entities, chains, residue source IDs | Missing writer |
| `struct_ref` | entity dbrefs and subchains | Missing writer |
| `chem_comp` | residue names and entity sequences | Missing writer |
| `exptl`, `diffrn`, `reflns`, `refine`, `tls`, `software` | typed metadata rows | Missing writer |
| `title_keywords` | typed metadata fields | Missing writer |
| `ncs`, `origx` | NCS and ORIGX fields | Missing writer |
| `struct_conf`, `struct_sheet` | helices and sheets | Missing writer |
| `struct_biol`, `assembly` | REMARK 300 detail and assemblies | Missing writer |
| `conn`, `cis`, `modres` | connections, cis-peptides, modified residues | Missing writer with reader gaps noted below |
| `atom_type` | atom elements | Missing writer |

`MmcifOutputGroups(true)` sets every group to true except `auth_all`, which is
false. The Rust default must reproduce that constructor, not derive all booleans
as true mechanically.

## Source-Model Gaps That Block Exact Writing

### Modified-residue identity is discarded

Gemmi stores `ModRes::res_id.name` and writes it to both
`auth_comp_id` and `label_comp_id`. `read_struct_mod_residue` currently reads
that value into `_mod_name` and discards it. `BioModRes` has no corresponding
field. This is a real roundtrip data-loss bug and must be fixed by modeling the
residue name, not by reconstructing it heuristically from hierarchy rows.

### Connection types are artificially narrowed

Gemmi's modeled connection type set includes `covale`, `disulf`, `hydrog`,
`metalc`, and `unknown`. `BioConnectionType` currently contains only the first,
second, and fourth; the reader rejects hydrogen bonds and unknown/null values.
The enum and source mappings must be completed before `write_struct_conn` can
be exact.

### Hybrid-36 input remains unported

Gemmi `read_pdb` accepts hybrid-36 atom serial and residue sequence fields.
`parse_decimal_i32` currently rejects alphabetic fields with an unsupported
error. This is a valid PDB source branch and remains a reader port gap.

### mmCIF identifiers are narrower than Gemmi strings

Gemmi chain, subchain, and entity identifiers are `std::string`. COSMolKit's
`PdbChainId([u8; 4], u8)` rejects mmCIF chain identifiers longer than four
bytes. This type currently conflates legacy-PDB fixed-width identifiers with
mmCIF identifiers. Exact complete-mmCIF behavior requires a source-identifier
type that preserves arbitrary valid mmCIF strings while retaining PDB field
validation at the PDB parser boundary. It must not truncate, hash, or invent an
alias.

### Category-presence metadata differs from typed optional fields

Gemmi metadata records whether a field was present independently from its
default value. COSMolKit generally uses `Option`, which preserves this
distinction, but `BioAnisotropicB` and `BioTlsGroup` use NaN sentinels. Writer
branches must inspect exactly the corresponding option/sentinel state and
must not emit optional columns merely because another row uses them unless
Gemmi does so.

## Existing Reader Markers And Misleading Guards

The remaining three `Gemmi❌❌` and ten `Gemmi❗✔️` markers are concentrated in
the flat-row PDB builder. Most describe implementation-shape differences, not
an observable unsupported branch. They must be re-audited after the writer and
reader gaps above; behaviorally covered lines should receive honest markers,
while genuinely different behavior must remain non-green.

`BioPdbReadParams::reject_unported_records` also contains stale rejection
branches for TER and connection records whose normal parser branches now
exist. Since its default is false, it does not affect ordinary reads, but it
misstates the current port boundary and creates a second behavior mode not
present in Gemmi. It should be removed or narrowed only after source comparison
and focused tests.

## Bio Operation Runtime

`BioOpParts::clear_cache()` currently records trace state but clears no storage.
This is not a hidden cache-preservation implementation: `BioStructure` has no
derived-cache block at all. The contract currently declares cache classes that
the data model does not store. The source writer does not require adding such
caches. The correct closure is to keep the operation contract aligned with
real storage, rather than invent empty cache values merely so `clear_cache()`
has something to mutate. This is an operation-system cleanup, not part of the
Gemmi writer call path.

## Required Completion Evidence

Completion requires all of the following:

1. one CIF document implementation supports both parsed data and writer
   mutation/serialization without a parallel ad-hoc text builder;
2. every source writer function is copied beside its Rust counterpart with
   line-by-line Gemmi two-axis markers;
3. every represented output group is exact and independently testable;
4. source-model gaps above are fixed at the typed model boundary;
5. public string and file APIs are added only after the full internal path;
6. pinned-Gemmi small-fixture output is compared exactly and reparsed
   semantically;
7. Python exposes only the completed `BioStructure` writer;
8. current policy, docs, examples, and the wwPDB experiment reflect the actual
   implemented boundary; and
9. the final audit lists any Gemmi field still not represented without
   claiming complete preservation for it.
