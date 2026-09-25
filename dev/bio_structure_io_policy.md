# BioStructure PDB/mmCIF IO Policy

This policy defines the boundary between Gemmi structural IO and RDKit molecule
compatibility. Ownership follows [crate_architecture.md](./crate_architecture.md);
execution and acceptance belong only to the
[split-crate plan](./plans/crate_architecture_completion_plan.md).
This document does not claim that a reader, writer, or conversion is implemented.

## Models and owners

| Model | Preserved state and role |
|---|---|
| `BioStructure` | Complete modeled structural hierarchy, mixed residue kinds, entities, coordinates, assemblies, and crystallographic metadata |
| `Protein` | Explicit amino-acid projection and ergonomic traversal; not lossless structural IO |
| `Molecule` | Chemical graph, chemical state, properties, and conformers; not the full biomolecular hierarchy |

Structural values and Protein projection belong to `cosmolkit-bio`.
Gemmi parsing, format handling, and structural serialization belong to
`cosmolkit-io`. RDKit molecule conversion also belongs to IO, reusing core
chemistry algorithms. Only `cosmolkit` constructs or accepts live molecules.

Read-only structural child views must not deep-clone the complete structure.
Structural mutation follows
[bio_structure_operation_contract_design.md](./bio_structure_operation_contract_design.md)
and [policy_invariants.md](./policy_invariants.md); it is not permission to
edit hierarchy rows through an unrestricted public storage interface.

## One structural reader

The data path is:

```text
PDB/mmCIF text -> Gemmi-aligned parsing -> BioStructure
                                          |
                                explicit conversion profile
                                          |
                             detached chemical graph + postprocessing
                                          |
                              cosmolkit validated Molecule construction
```

There is no competing public RDKit structural parser. Do not expose parallel
Gemmi/RDKit parser modules or place RDKit molecule rules inside the Gemmi parser.
RDKit's PDB reader specifies molecule compatibility behavior; it is not a
replacement structural reader or a pure subset of Gemmi.

Gemmi sources include `src/pdb.cpp`, `src/mmcif.cpp`,
`include/gemmi/mmcif.hpp`, and `include/gemmi/mmread.hpp` under the pinned
`third_party/gemmi` checkout. The structural scope includes:

- model/chain/residue/atom order and hierarchy; serials, names, insertion codes;
- altloc, occupancy, B factors, formal charge, XYZ and ANISOU;
- header records, entity/sequence/DBREF and author/label source identities;
- HELIX, SHEET, SSBOND, LINK, CISPEP, MODRES, TER and connection records;
- CRYST1, SCALE, ORIGX, MTRIX, cell, space group, NCS and assembly metadata;
- mmCIF atom-site, entity, sequence, connectivity and crystallographic categories.

A complete declared source profile must be implemented and validated before
acceptance. This inventory is not permission to skip a required category.

## Public construction and identifiers

Public naming, defaults and format/error dispatch are declared in
[`crates/cosmolkit/src/binding_contract/registry.rs`](../crates/cosmolkit/src/binding_contract/registry.rs)
under the [public API design](./public_api_design.md).
The split-crate structural text constructors are model-qualified functions in
`cosmolkit::bio`, not a second public parser and not inherent IO methods added
to an externally owned `BioStructure`. Neither a `bio -> io` reverse
dependency nor an illegal cross-crate inherent implementation is acceptable.

The old constructor spellings and old core IO paths are historical locators,
not canonical aliases. Registration precedes public exposure. In particular,
a subset-reader name must not be used to imply complete structural reading.

The declared source-ID boundary uses at-most-four-byte `PdbChainId` values
for chain/subchain hierarchy references. Over-width identifiers return the
documented structured boundary error, never truncation, hashing or invented
aliases. Entity source identifiers remain variable-length strings. Author and
label chain identities remain distinct. Exact AtomName whitespace and length,
altloc, source order and coordinate behavior follow the frozen source/profile
mapping; no convenient lossy normalization is allowed.

## Molecule conversion

Conversion is explicit because structural and chemical models preserve
different invariants. Its contract must identify:

- retained/filtered atoms and source model/altloc selection;
- preserved, synthesized and intentionally lost fields;
- source identifiers and PDB residue information attached to molecule atoms;
- imported CONECT/mmCIF connections versus source-backed proximity bonding;
- atom/residue mappings, sanitization and hydrogen-removal order;
- structured errors for independently unmodeled capabilities.

Do not guess chemistry or relabel supported-input failures as unsupported.
A protein-only projection cannot claim `Chem.MolFromPDBBlock()` equivalence:
it can discard ligands, water, ions, nucleic acids and connection context.

Detailed filtering and atom/bond behavior must be reproduced from the pinned
RDKit PDB reader and its reached helpers, with source anchors in the owning
implementation and regression coverage for the conversion profile above.
The facade may expose destination/source-oriented conversion APIs only through
the canonical public contract. Historical `RdkitPdbMolProfile`,
`PdbMoleculeConversionError`, old `_with_options` examples and free-function
shims do not authorize compatibility aliases or override current naming rules.

There is no independent RDKit mmCIF parser. Molecule input from mmCIF applies
the documented conversion profile after the Gemmi structural parse.

## Structural writers

Structural serialization targets `BioStructure` and uses one Gemmi-aligned
document builder and CIF serializer. A PDB-to-mmCIF workflow is composition of
the public reader and writer, not a separate format-pair parser.

A writer emits represented structural state. It cannot promise preservation of
unmodeled/private CIF categories unless the model explicitly retains them.
Canonical writer names, parameter types and language projections must be
registered under [public_api_design.md](./public_api_design.md), not copied
from historical signatures.

Do not expose `Protein::to_mmcif()` as lossless structural conversion or
`Molecule::to_mmcif()` as complete structural serialization. RDKit-based
`Molecule` PDB writing is a separate chemical-graph writer, not structural IO.

## Completion and change boundaries

The structural reader must close its declared Gemmi profile. Molecule input
must additionally close the conversion and RDKit postprocessing profile.
These are separate acceptance obligations, not independent execution queues.

Source scope, test evidence and unresolved dependencies belong to the owning
unit reports. The sole plan schedules structural prerequisites before their
consumers. A task that genuinely requires a parallel parser, new ownership
boundary, or changed model support limit must stop for explicit approval.
