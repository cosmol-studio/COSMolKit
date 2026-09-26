# Coordinate Storage and Selection Contract

Status: approved design. Approval date: 2026-09-26.

## Storage and purpose

CK stores 2D layouts and 3D conformers in separate typed collections.
They may coexist on one molecule. A 2D layout supports depiction; a 3D
conformer supports spatial calculations. Neither collection is a fallback for
the other. Coordinate queries are read-only and never generate coordinates.

There is no cross-dimension primary conformer or semantic insertion order.
A coordinate set is selected by dimension and its existing conformer ID,
not by vector position or the smallest ID. IDs are scoped to their dimension.
The same numeric ID in the two collections is not ambiguous once a dimension
is specified. Existing model validation rejects duplicate IDs within a dimension.

Coordinate generation changes its declared dimension and preserves the other.
Topology edits must map every retained coordinate set through the authoritative
atom mapping. These rules do not authorize a new layout-management API or a
wholesale change to existing conformer identity/remapping behavior.

Original file ordering is I/O provenance only when a concrete round-trip
requirement calls for it; it must not determine algorithm defaults. No such
provenance storage is required by this contract.

## Selection

A consumer declares which dimensions it accepts. Dimension-specific consumers
must not silently fall back to another dimension. For an exporter accepting
both dimensions, automatic selection means:

| Available coordinate sets | Automatic selection |
|---|---|
| None | Emit no coordinate block |
| Exactly one, in either collection | Use that set |
| More than one, including multiple sets of one dimension | Return a structured ambiguity error |

An explicit selection consists of dimension and ID. A missing selected ID
returns a structured error, even if some other coordinate set exists. Selection
does not generate, project, convert, reorder, or delete coordinates. It does
not reinterpret XYZ geometry as 2D merely because z is zero.

Coordinate-disabled export does not resolve selection or raise selection errors.
Existing structural validation still applies. All coordinate-dependent stages
of an export use the same selected set; text output and wedge geometry must
not independently choose different conformers.

Drawing uses an explicitly selected 2D layout; spatial algorithms use an
explicitly selected 3D conformer. Their convenience defaults must document a
dimension-specific choice rather than inherit a mixed collection's first entry.

## Approved reference difference

CK-COORD-002: CX export intentionally replaces RDKit's implicit first-conformer
selection with unique-or-explicit selection. This removes insertion-order and
dimension-preference ambiguity while preserving CK's separate 2D/3D model.

Pinned reference behavior is in RDKit CXSmilesOps.cpp get_coords_block /
getCXExtensions and ROMol.cpp getConformer(-1). Retain the verbatim source
anchors and mark the selection difference explicitly in the implementing
function; do not describe it as exact parity or an unsupported chemistry case.
Swapping mixed-dimension input order must not silently change CK's default
selection: both inputs are ambiguous unless a coordinate set is selected.

After selection, numerical formatting, atom-order mapping and coordinate-based
stereo/wedge calculations still follow the pinned source for that selected set.
This approval does not permit approximation in those algorithms or changes to
input interpretation, sanitization, stereo perception, or is_3d handling.

## CX writer implementation contract

The existing detached CxSmilesWriteParams gains coordinate_selection with
default CxCoordinateSelection::Auto. The domain enum is:

```rust
enum CxCoordinateSelection {
    Auto,
    TwoD { id: usize },
    ThreeD { id: usize },
}
```

The existing writer error channel carries typed ambiguity counts and typed
missing-selection data. Selection is resolved once for the coordinate-enabled
export and borrowed by both rendering and coordinate-dependent wedge processing.
No global ordering field, parser ID workaround or public mutable coordinate
storage is added.

This is a domain-crate API change, not a new runtime/binding export. Any future
runtime projection must first register its selector, defaults, errors and this
approved difference in binding_contract. No registry support state is promoted
by this design approval.

## Validation and scope

Use focused fixed regressions for zero/one/multiple sets, mixed input order,
same numeric ID across dimensions, noncontiguous IDs versus positions, missing
IDs, disabled coordinate fields and shared selection for rendering/wedges.
Preserve input values and earlier counterexamples. Compare selected-coordinate
algorithm behavior with the corresponding pinned reference conformer, not the
reference's implicit mixed-order choice.

This contract is authoritative for storage separation and selection. Other
proposals in coordinate_model_improvement_draft.md remain deferred. The sole
execution ledger and its delegated SMILES plan record implementation progress;
this document does not claim that every consumer already implements the design.
