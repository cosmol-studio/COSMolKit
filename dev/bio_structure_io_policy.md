# BioStructure PDB/mmCIF IO Policy

This document fixes the project-level boundary between Gemmi-derived structural
IO and RDKit-derived molecule compatibility behavior. The detailed execution
plan is [`pdb_mmcif_gemmi_primary_plan.md`](pdb_mmcif_gemmi_primary_plan.md).

The short rule is:

```text
BioStructure owns complete structural IO. Protein is an explicit amino-acid
projection. Molecule owns chemical-graph compatibility behavior.
```

COSMolKit must not expose two parallel public parser families for the same user
task. In particular, do not create competing public APIs such as
`pdb_parser_gemmi`, `pdb_parser_rdkit`, `mmcif_parser_gemmi`, and
`mmcif_parser_rdkit`.

## Public Model

The public structural model is COSMolKit's own `BioStructure` hierarchy:

- `BioStructure`
- `ModelRow`
- `ChainRow`
- `ResidueRow`
- `AtomRow`
- `CoordinateBlock`
- registered `BioStructure` operations

Public structural IO should read into this model. Mutation of this model must
follow the BioStructure operation rules in `dev/README.md`,
`dev/policy_invariants.md`, and `dev/topology_operations.md`.

The public object boundary is:

| Object | Preserved state | Owns | Must not imply |
|---|---|---|---|
| `BioStructure` | complete modeled hierarchy, entities, mixed residue kinds, coordinates, assemblies, and crystallographic metadata | structural PDB/mmCIF/mmJSON reading and future structural writing | RDKit molecule-graph semantics |
| `Protein` | an amino-acid-only projection of `BioStructure` | ergonomic protein chain/residue/atom traversal | lossless structural format conversion |
| `Molecule` | chemical graph, chemical state, properties, and conformers | RDKit-compatible molecule parsing, conversion, and molecule writers | preservation of the complete biomolecular hierarchy |

Rust exposes `BioStructure` directly. Python must expose the same complete
structural concept rather than forcing full-structure workflows through
`Protein`. Child model/chain/residue/atom/entity values are read-only shared
views; producing them must not deep-clone the complete structure.

## Canonical Structural IO Source

Gemmi is the canonical upstream source for PDB/mmCIF structure-file IO.

The Gemmi-aligned path owns:

- PDB/mmCIF coordinate hierarchy parsing
- model, chain, residue, atom, altloc, occupancy, B-factor, formal charge, and
  coordinate rows
- mmCIF `_atom_site` column handling
- structural unsupported-feature boundaries
- entity, assembly, secondary-structure, symmetry, and crystallographic
  metadata represented by `BioStructure`

The current Rust implementation of this path lives in:

```text
crates/cosmolkit-core/src/io/bio.rs
```

The canonical Rust structural read methods are:

```text
BioStructure::from_pdb
BioStructure::from_mmcif
BioStructure::from_structure_str
BioStructure::from_str_with_format
BioStructure::from_pdb_str
BioStructure::from_pdb_str_with_params
BioStructure::from_mmcif_str
```

`BioStructure::from_structure_str(...)`,
`BioStructure::from_str_with_format(...)`, and
`BioStructure::from_pdb_str(...)` are the canonical public dispatch entry
points for Gemmi-aligned structural reads. Format helpers are exposed as
`BioStructure` constructors so callers choose the target data model explicitly
instead of invoking ambiguous free functions.

`read_mmcif_atom_site_subset_from_str(...)` is a deprecated historical name.
The implementation now reads the declared Gemmi-aligned structural surface,
not only `_atom_site`, so it must not be presented as the primary API.

## RDKit PDB/mmCIF Scope

RDKit PDB behavior is not the canonical structural reader for COSMolKit, and it
must not become an independent text parser that competes with the Gemmi path.

RDKit-derived code may be ported only for Molecule compatibility tasks after
the PDB/mmCIF text has gone through the Gemmi-primary structural path, for
example:

- RDKit-compatible PDB block writing from `Molecule`
- RDKit-compatible `Molecule::from_pdb_block` behavior built on
  `BioStructure` conversion
- RDKit-compatible `AtomPDBResidueInfo` attachment to molecule atoms
- RDKit-compatible bond perception, sanitization, and residue metadata behavior
  for molecule graph construction

RDKit-derived PDB code must output or operate on `Molecule`, not a competing
protein hierarchy. The molecule readers are layered over the Gemmi-primary
`BioStructure` parse rather than a separate public structural parser.

Current primary API shape:

```text
Molecule::from_pdb_block(input)
Molecule::from_pdb_block_with_options(input, StructureMoleculeOptions)
Molecule::from_mmcif_block(input)
Molecule::from_mmcif_block_with_options(input, StructureMoleculeOptions)
BioStructure::to_molecule()
BioStructure::to_molecule_with_options(StructureMoleculeOptions)
Molecule::to_pdb_block(...)
```

`StructureMoleculeOptions` and `StructureMoleculeConversionError` are named for
the model boundary because the same conversion policy applies after PDB or
mmCIF structural parsing. Historical `RdkitPdbMolProfile` and
`PdbMoleculeConversionError` names are deprecated compatibility aliases.
Top-level conversion functions are likewise compatibility shims; new code
should use methods on the destination or source value.

Unacceptable future API shape:

```text
pdb_parser_rdkit::read_structure(input)
mmcif_parser_rdkit::read_biostructure(input)
```

There is no independent RDKit mmCIF parser. `Molecule::from_mmcif_block()`
starts from the Gemmi-primary structural parse and deliberately applies the
documented structure-to-molecule conversion options.

## Conversion Boundary

Conversions between `BioStructure` and `Molecule` must be explicit because the
two models preserve different invariants.

`BioStructure` preserves structural hierarchy and source identifiers. `Molecule`
preserves chemical graph state, valence/sanitization state, conformers, and
atom/bond metadata.

The conversion API must define:

- which structural fields are preserved
- which molecule fields are synthesized
- whether bonds are imported, inferred, or unsupported
- whether sanitization runs
- how atom/residue mappings are reported
- what information is intentionally lost

Do not hide conversion behind structural parser names. A PDB/mmCIF structural
read is not the same operation as molecule graph construction.
`Molecule::from_pdb_block` is allowed only as a molecule compatibility API
whose internal pipeline is Gemmi structural parse, explicit conversion, then
RDKit-compatible molecule post-processing.

## Structural Writer Boundary

Structural writers belong only to `BioStructure`. A PDB-to-mmCIF workflow is
composition through the public model, not a separate format-pair function:

```text
BioStructure::from_pdb_str(input)?.to_mmcif()
```

The following API is planned but is not implemented yet:

```text
BioStructure::to_mmcif()
BioStructure::to_mmcif_with_options(MmcifWriteOptions)
BioStructure::write_mmcif(path)
```

Do not add placeholders that return unsupported errors merely to make these
methods appear present. They become public only after Gemmi
`make_mmcif_document`, `update_mmcif_block`, output-group behavior, and CIF
serialization have been source-ported and tested.

Do not add `Protein::to_mmcif()`: `Protein` has already discarded non-protein
rows, so that name would look like lossless format conversion while silently
writing a projection. Do not add `Molecule::to_mmcif()`: `Molecule` does not
represent the complete structural hierarchy or mmCIF metadata.

The initial structural writer will generate a canonical mmCIF document from
the state modeled by `BioStructure`. It must not claim to preserve arbitrary
unmodeled or private categories from an input CIF document unless the public
model later retains that source document explicitly.

## Planning Rule

When porting README-claimed PDB/mmCIF/protein behavior:

1. Use `io::bio` and `BioStructure` for all PDB/mmCIF reading.
2. Use RDKit source only for Molecule compatibility behavior layered after
   `BioStructure` exists.
3. Keep unsupported branches explicit.
4. Do not expose quarantined `pdb_parser.rs` or `mmcif_parser.rs` as public
   modules.
5. Do not merge Gemmi and RDKit code into one large mixed parser. Keep the
   Gemmi parser source-scoped and keep RDKit rules in the conversion or molecule
   compatibility layer.

This policy is binding for future agents. If a task appears to require a public
parallel parser API, stop and ask the human author before implementing it.
