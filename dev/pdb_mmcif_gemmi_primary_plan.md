# PDB/mmCIF Gemmi-Primary Porting Plan

This plan replaces any idea of maintaining separate Gemmi and RDKit PDB/mmCIF
readers.

The target architecture is:

```text
PDB/mmCIF text
  -> Gemmi-aligned structural parser
  -> BioStructure
  -> explicit BioStructure-to-Molecule conversion profile
  -> RDKit-compatible molecule graph post-processing
  -> Molecule
```

There must be no second public RDKit parser path that reads PDB/mmCIF text
directly into a molecule graph independently of the Gemmi-aligned structural
reader.

## Design Decision

PDB and mmCIF reading are Gemmi-primary.

RDKit's PDB reader is treated as a compatibility specification for molecule
graph construction, not as a competing structural parser. The project ports
Gemmi's structural behavior first, then overlays the RDKit-specific molecule
construction behaviors needed by `Molecule` APIs.

This is more precise than saying RDKit PDB reading is a pure Gemmi subset. The
structure fields mostly belong to the Gemmi side, but RDKit has molecule-graph
behavior that Gemmi does not own:

- alternate location and dummy atom filtering flags
- multi-model handling flags
- `AtomPDBResidueInfo` attachment
- CONECT-to-bond rules
- proximity bonding
- sanitize and remove-H behavior
- bond order repair and zero-bond fallbacks

These behaviors belong after structural parsing, during explicit conversion and
RDKit-compatible molecule post-processing.

## Layer 1: Gemmi-Aligned Structural Parse

Goal: build a complete `BioStructure` from PDB/mmCIF input.

Source target:

```text
third_party/gemmi/src/pdb.cpp
third_party/gemmi/src/mmcif.cpp
third_party/gemmi/include/gemmi/mmcif.hpp
third_party/gemmi/include/gemmi/mmread.hpp
```

This layer owns:

- model, chain, residue, atom hierarchy
- atom serial, atom name, residue name, residue id, insertion code, chain id
- altloc, occupancy, B factor, formal charge, coordinates, ANISOU
- PDB header and metadata records when modeled
- entity, sequence, DBREF, struct asym, polymer metadata
- HELIX, SHEET, SSBOND, LINK, CISPEP, MODRES and related structural records
- CRYST1, SCALE, ORIGX, MTRIX, cell, space group and NCS metadata
- mmCIF `_atom_site` plus entity, sequence, connectivity, assembly and
  crystallographic categories

Current state: `crates/cosmolkit-core/src/io/bio.rs` implements the coordinate
hierarchy plus PDB `SEQRES`, selected `HEADER`/`TITLE`/`KEYWDS`/`EXPDTA`/
`AUTHOR` metadata, `CRYST1`, and mmCIF `_entity`, `_entity_poly`,
`_entity_poly_seq`, and `_struct_asym` entity/sequence binding. This is still
not complete Gemmi structural parse: DBREF/struct_ref, TER splitting, HELIX,
SHEET, SSBOND, LINK, CISPEP, MODRES, CONECT, assemblies, cell transforms beyond
CRYST1, and broad mmCIF category handling remain unported.

## Layer 2: COSMolKit Conversion Profile

Goal: convert `BioStructure` to `Molecule` explicitly, with no hidden parser
semantics.

This layer is COSMolKit-owned. It must define:

- which atoms become molecule atoms
- how residue and chain source ids map to `AtomPDBResidueInfo`
- whether bonds come from PDB CONECT, mmCIF connectivity, proximity perception,
  or an explicit unsupported error
- how multiple models are selected or split
- how altlocs are selected or rejected
- how atom/residue mappings are reported
- what information is intentionally lost when moving from `BioStructure` to
  `Molecule`

No conversion API may silently guess chemistry. Unsupported conversion branches
must return structured errors.

## Layer 3: RDKit-Compatible Molecule Behavior

Goal: implement `Molecule` APIs that match RDKit where COSMolKit chooses to
provide RDKit compatibility.

Source target:

```text
third_party/rdkit/Code/GraphMol/FileParsers/PDBParser.cpp
third_party/rdkit/Code/GraphMol/FileParsers/PDBWriter.cpp
```

This layer may define APIs such as:

```text
Molecule::from_pdb_block_with_params(input, params)
mol_to_pdb_block(molecule)
```

The read API must use the Gemmi-primary structural parse layer instead of
introducing a separate parser. RDKit source is used to reproduce molecule graph
decisions after `BioStructure` exists.

PDB writing from `Molecule` may remain RDKit-primary because it is a
`Molecule -> PDB block` operation, not structural file reading. If COSMolKit
later adds `BioStructure -> PDB/mmCIF` writing, that writer should be
Gemmi-primary.

## Forbidden Designs

Do not implement or expose:

```text
pdb_parser_rdkit::read_structure(input)
pdb_parser_rdkit::read_molecule_independent_of_biostructure(input)
mmcif_parser_rdkit::read_biostructure(input)
```

Do not place RDKit molecule construction rules inside the Gemmi structural
parser. Keep the parser source-mapped to Gemmi and keep RDKit behavior in the
conversion or molecule compatibility layer.

Do not expose `crates/cosmolkit-core/src/unported/pdb_parser.rs` or
`crates/cosmolkit-core/src/unported/mmcif_parser.rs` as public modules. They remain
quarantine placeholders unless they are repurposed into private compatibility
helpers that still obey this plan.

## Execution Order

1. Expand `io::bio` toward complete Gemmi-aligned PDB/mmCIF `BioStructure`
   parsing.
2. Add typed `BioStructure` fields only when required by source-backed Gemmi
   behavior or by an explicit conversion profile.
3. Define a `BioStructure -> Molecule` conversion profile with explicit errors.
4. Implement RDKit-compatible `Molecule::from_pdb_block` behavior as a thin API
   over the Gemmi-primary parse plus conversion plus RDKit post-processing.
5. Keep RDKit PDB writer work under the existing `Molecule -> PDB block` path.

Completion of PDB/mmCIF reading means the Gemmi structural layer is complete for
the declared support profile. Completion of RDKit PDB molecule input means the
Gemmi layer and conversion layer satisfy the declared RDKit compatibility
profile.
