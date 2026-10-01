# BIO public objects and lightweight operations

`cosmolkit` owns the public `BioStructure` and `Protein` objects. Both expose
associated PDB/mmCIF text constructors. `Protein` is an explicit amino-acid
projection, not a lossless substitute for a mixed structural hierarchy.

`cosmolkit-bio` owns detached structure data, row/identifier types, local
validation and structural algorithms. `cosmolkit-io` owns the sole
Gemmi-aligned structural readers and writers; it returns detached data, never
public objects. Neither lower crate depends on `cosmolkit`. Public adapters
only call these owners and wrap validated results. Molecule conversion remains
a separate explicitly named IO algorithm plus the molecule runtime commit.

The public `bio` feature enables `cosmolkit-io/bio`. IO's BIO dependency,
structural PDB/mmCIF readers and CIF parser are optional behind that feature;
ordinary `io` and `serialization` do not enable them. The default public `full`
bundle includes `bio`. BIO-specific IO tests explicitly require `bio`.

## Data and sharing

The detached structure stores independent shared blocks for models, chains,
residues, atoms, entities, coordinates, connections, cis-peptides, modified
residues, helices, sheets, metadata, source state, crystal data, NCS operators
and assemblies. Input format is a copied value. There is one representation,
not a second chemistry model in the public crate. Local structure validation
checks hierarchy spans, references and coordinate alignment. Detached data may
be edited by algorithms; public objects never expose mutable data or DerefMut.

Cloning a public object shares blocks. A generated operation creates a working
snapshot and performs COW on its declared writable blocks. Read blocks are
borrowed. Unaffected blocks remain shared. Retaining the original snapshot for
rollback means writable blocks can require copying even for an originally
unique in-place receiver; no zero-copy in-place guarantee is made.

## One declaration and one implementation

The lightweight `bio_structure_ops!` declaration generates operations for both
public types from the same mechanism. An entry names its target types, value
method, trailing-underscore in-place method, implementation, arguments and
read/write fields. It generates the access struct, both methods and metadata.
The implementation receives only declared `&Block` / `&mut Block` references,
not the public object, complete detached data or commit authority. Algorithms
remain below the public crate. There is no parallel handwritten permission list.

Both methods run the same implementation on a working snapshot. Final
structural validation precedes return or replacement of the receiver. Protein
validation additionally preserves the amino-acid-only invariant. Errors and
unwinding panics leave the receiver and shared snapshots unchanged; aborting
panics terminate the process. There is no tracing of read/write use, mandatory
write-back accounting, molecule-derived cache protocol or runtime permission
checker. Declaring write permission does not require using it. Compile-time
field availability supplies the access boundary.

Text construction has no in-place partner and belongs only in the public
binding registry, not the mutation registry. Public callables and supporting
types have compiler-checked binding entries. New entries remain experimental
until explicitly reviewed; passing tests does not promote status. Python/JS
names describe future projections, not implemented adapters.

## Source and format contracts

Selection copying follows pinned Gemmi `Selection::copy_selection` and
`Structure::empty_copy`. Retain source-assigned assembly and other metadata,
including chain/subchain names absent from the selected hierarchy. These names
are source metadata, not local row references; construction must not reject,
prune or rewrite them merely because selection omitted rows. A consumer needing
actual targets must resolve names and explicitly handle absence. Local row IDs,
parent/spans, entity-row references and coordinate alignment remain validated.
Accepted parents remain present even when no children match. Fields not assigned
by `empty_copy` follow its defaults; this is not whole-object cloning or the
separate Protein projection contract. This rule does not authorize assembly
expansion, symmetry calculations or an alternate validation-bypass constructor.

Preserve the pinned reader's order, options, errors, source identifiers,
metadata, coordinates, altlocs, crystal/NCS and assemblies across its declared
scope. Existing representational limits return structured errors, never
padding, truncation, aliases or guessed values. Author and label identities
remain distinct. Do not expose competing RDKit and Gemmi structural readers.
RDKit molecule compatibility follows explicit conversion after structural IO.
Structural writers cannot claim preservation of unmodeled categories; Protein
serialization must not claim lossless structural conversion. Source/performance
debts remain open independently of public API plumbing.

## Validation

Fixed local regressions cover both objects' constructor/projection behavior,
registry signatures, options and error sources; operation tests cover matching
value/in-place results, unchanged snapshots/block sharing and rollback after
body or final-validation failure. Compile-fail cases use real operation modules
to reject undeclared fields and direct object storage access. No corpus runner
is introduced into domain crates. Corpus comparison belongs in parity-tests.
