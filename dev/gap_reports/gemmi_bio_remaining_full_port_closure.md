# Gemmi BioStructure Remaining Full-Port Closure

## Audit Boundary

This closure audit compares the completed `BioStructure` reader and mmCIF
writer line with Gemmi commit
`5cc1c23c6007e0e6cbd69289c6f7c0bff50e943e`.

The source boundary reviewed here is:

- `third_party/gemmi/src/pdb.cpp`
- `third_party/gemmi/include/gemmi/pdb.hpp`
- `third_party/gemmi/src/mmcif.cpp`
- `third_party/gemmi/include/gemmi/mmcif.hpp`
- `third_party/gemmi/include/gemmi/mmcif_impl.hpp`
- `third_party/gemmi/src/to_mmcif.cpp`
- `third_party/gemmi/include/gemmi/to_mmcif.hpp`
- `third_party/gemmi/include/gemmi/to_cif.hpp`
- `third_party/gemmi/include/gemmi/cifdoc.hpp`

The public ownership boundary remains unchanged. Structural PDB/mmCIF/mmJSON
reading and mmCIF writing belong to `BioStructure`. `Protein` remains an
amino-acid projection and `Molecule` remains a chemical graph; neither owns a
parallel structural writer.

## Closed Source Paths

The implemented line now contains one source-backed path for each represented
operation:

- PDB record parsing, finalization, hybrid-36 serial and sequence decoding,
  metadata, crystallographic transforms, entities, connections, secondary
  structure, assemblies, and anisotropic displacement values;
- mmCIF and mmJSON document parsing, structure construction, atom/entity
  population, metadata, connections, assemblies, secondary structure,
  chemical-component dispatch, and format detection;
- mmCIF document mutation, category replacement, loop construction, CIF
  quoting and serialization;
- all 34 `MmcifOutputGroups` switches with Gemmi defaults;
- `add_cif_atoms`, assemblies, NCS operations, structural connections,
  cis-peptides, represented metadata categories, and the public Rust/Python
  `BioStructure` writer entry points.

The former reader gaps from the opening inventory are closed: modified-residue
identity is typed, all Gemmi connection kinds are represented, PDB hybrid-36
fields are decoded through the Gemmi algorithm, stale record-rejection options
are removed, and the flat-row builder uses indexed per-chain residue lookup.

## Marker Audit

The implementation files contain no non-green Gemmi source marker:

| File | `Gemmi✔️✔️` occurrences | Non-green implementation markers |
|---|---:|---:|
| `crates/cosmolkit-core/src/io/bio.rs` | 1,397 | 0 |
| `crates/cosmolkit-core/src/io/bio/cif.rs` | 146 | 0 |
| `crates/cosmolkit-core/src/io/bio/sprintf.rs` | 404 | 0 |
| `crates/cosmolkit-core/src/io/bio/writer.rs` | 1,002 | 0 |

The strings `Gemmi❌❌` and `Gemmi❗✔️` at the top of `io/bio.rs` define the
marker notation. They are not attached to a source line or implementation and
therefore are not remaining port markers.

## Structured Rejection Audit

Five calls to the local `unsupported(...)` constructor remain in the reader.
Four describe malformed or inapplicable input rather than an unported Gemmi
branch:

| Reader condition | Classification |
|---|---|
| no `_atom_site` data in a structure block | required structural input is absent |
| an empty mmCIF document | no first block exists |
| a chemical-component block has no recognized coordinate triplet | no selectable Gemmi chemical-component model exists |
| a chemical-component block has no `_chem_comp_atom.atom_id` loop | required chemical-component input is absent |

These exits do not return a plausible partial structure and do not conceal an
alternative implementation. The fifth exit is the explicit modeled-state
boundary described below.

## Remaining Unmodeled Source State

The completed line reproduces Gemmi behavior for state representable by
`BioStructure`. It does not claim that `BioStructure` is a byte-preserving
container for every possible CIF value or category.

### Arbitrary-length mmCIF hierarchy identifiers

Gemmi stores chain and subchain identifiers as `std::string`.
`BioStructure` currently stores author-chain and label-asym/subchain values in
`PdbChainId([u8; 4], u8)`. An mmCIF chain or subchain identifier longer than
four bytes therefore returns a structured unsupported error. The reader does
not truncate, hash, or invent an alias. Entity source identifiers themselves
are retained as `String`; the limitation is the hierarchy references carried
through `PdbChainId`.

Closing this boundary requires a dedicated variable-length structural
identifier type and a model-wide migration. It is not an empty parser branch
and must not be hidden behind a local conversion rule.

### Fixed-width atom and residue names

`AtomName` stores four bytes and `ResidueName` stores four bytes plus length,
matching the PDB-oriented public row model. Gemmi uses `std::string` for these
names. Consequently, atom or component names wider than the current row types
are not completely representable. Complete support requires variable-length
typed names across the structural model, lookup code, Python views, and writer;
it must not be implemented by a parser-only alias or hash.

### CIF categories outside the typed structural model

The parser consumes the source categories needed to construct the modeled
hierarchy and metadata. `BioStructure` does not retain the complete original
`cif::Document`, arbitrary private categories, unknown columns, comments, or
the original pair/loop layout. Its writer constructs canonical mmCIF from
typed state. Exact preservation of such document syntax or unmodeled category
payload is therefore outside the writer boundary.

### Source state intentionally represented at a different layer

`BioStructure` does not acquire RDKit molecule-graph state such as perceived
bonds, valence caches, sanitization state, or conformers. Those belong to the
explicit `BioStructure` to `Molecule` conversion boundary. Likewise,
`Protein` intentionally discards non-amino-acid rows and is not a lossless
structural serialization model.

## Closure Result

There are no remaining placeholder bodies, dispatch-only branches,
`unimplemented!` calls, or non-green Gemmi implementation markers in the
audited Bio reader/writer line. The implemented Gemmi call paths are complete
for the current typed `BioStructure` state. The remaining limitations above
are data-model boundaries, not silent alternate algorithms, and must be
addressed by explicit type evolution before their public support boundary can
expand.
