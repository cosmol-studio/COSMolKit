# BioStructure PDB/mmCIF IO Policy

This document fixes the project-level boundary between Gemmi-derived structural
IO and RDKit-derived molecule compatibility behavior. The detailed execution
plan is [`pdb_mmcif_gemmi_primary_plan.md`](pdb_mmcif_gemmi_primary_plan.md).

The short rule is:

```text
Gemmi-primary reading, one public structural model, RDKit-compatible molecule
post-processing.
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

## Canonical Structural IO Source

Gemmi is the canonical upstream source for PDB/mmCIF structure-file reading.

The Gemmi-aligned path owns:

- PDB/mmCIF coordinate hierarchy parsing
- model, chain, residue, atom, altloc, occupancy, B-factor, formal charge, and
  coordinate rows
- mmCIF `_atom_site` column handling
- structural unsupported-feature boundaries
- future structure-specific records such as entity, assembly, secondary
  structure, symmetry, and crystallographic metadata

The current Rust implementation of this path lives in:

```text
crates/cosmolkit-core/src/io/bio.rs
```

The current public functions are:

```text
read_pdb_coordinate_subset_from_str
read_pdb_coordinate_subset_from_str_with_params
read_mmcif_atom_site_subset_from_str
```

These names are intentionally subset-scoped. They must not be renamed into
complete PDB/mmCIF support until the code actually implements complete support.

## RDKit PDB/mmCIF Scope

RDKit PDB behavior is not the canonical structural reader for COSMolKit, and it
must not become an independent text parser that competes with the Gemmi path.

RDKit-derived code may be ported only for Molecule compatibility tasks after
the PDB/mmCIF text has gone through the Gemmi-primary structural path, for
example:

- RDKit-compatible PDB block writing from `Molecule`
- future RDKit-compatible `Molecule::from_pdb_block` behavior built on
  `BioStructure` conversion
- RDKit-compatible `AtomPDBResidueInfo` attachment to molecule atoms
- RDKit-compatible bond perception, sanitization, and residue metadata behavior
  for molecule graph construction

RDKit-derived PDB code must output or operate on `Molecule`, not a competing
protein hierarchy. If an RDKit-compatible PDB molecule reader is added later,
its public name must communicate molecule graph intent, and its implementation
must be layered over the Gemmi-primary `BioStructure` parse rather than a
separate parser.

Acceptable future API shape:

```text
Molecule::from_pdb_block_with_params(input, params)
mol_to_pdb_block(molecule)
```

Unacceptable future API shape:

```text
pdb_parser_rdkit::read_structure(input)
mmcif_parser_rdkit::read_biostructure(input)
```

RDKit mmCIF parser work is deferred unless a concrete Molecule compatibility
need is identified. Even then, the read path must start from the Gemmi-primary
structural parse unless the human author approves a design exception.

## Conversion Boundary

Conversions between `BioStructure` and `Molecule` must be explicit because the
two models preserve different invariants.

`BioStructure` preserves structural hierarchy and source identifiers. `Molecule`
preserves chemical graph state, valence/sanitization state, conformers, and
atom/bond metadata.

Any future conversion API must define:

- which structural fields are preserved
- which molecule fields are synthesized
- whether bonds are imported, inferred, or unsupported
- whether sanitization runs
- how atom/residue mappings are reported
- what information is intentionally lost

Do not hide conversion behind parser names. A PDB/mmCIF structural read is not
the same operation as molecule graph construction. `Molecule::from_pdb_block`
is allowed only as a molecule compatibility API whose internal pipeline is
Gemmi structural parse, explicit conversion, then RDKit-compatible molecule
post-processing.

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
