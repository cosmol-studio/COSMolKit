# Coordinate Model Improvement Draft

Status: deferred proposal, not an approved architecture or execution plan.

Date: 2026-09-16.

Revisit this document at the final review of the current migration work. Do
not interrupt, reorder, or expand the active migration based on this draft.
`dev/crate_architecture.md` remains the architectural authority and
`dev/plans/crate_architecture_completion_plan.md` remains the sole execution
queue. Any eventual semantic changes require explicit approval and integration
into the authoritative documents before implementation.

## 1. Assessment

Separating 2D depictions from 3D conformers addresses a real usability problem.
RDKit stores both in its conformer collection, and `Compute2DCoords` defaults
to clearing existing conformers. Retaining them is possible, but callers must
manage mixed dimensions and select the intended conformer. See the official
[C++ introduction](https://www.rdkit.org/docs/GettingStartedInC%2B%2B.html)
and [depiction API](https://www.rdkit.org/docs/source/rdkit.Chem.rdDepictor.html).

COSMolKit's dimension-specific accessors are a useful improvement:

```rust
mol.coordinates_2d();
mol.conformers_3d();
```

Separate storage can prevent depiction generation from destroying spatial
conformers and avoid ambiguous coordinate queries. Typed XY/XYZ rows, read-only
access, explicit value-returning transformations, and runtime-controlled atom
mapping are worth retaining.

However, separate storage alone does not establish preservation, selection,
or stereochemical guarantees. At the review date, the main depiction and
conformer generation entrypoints were still unsupported placeholders. This
draft describes desired contracts, not verified delivery of those features.

The main problem is overlapping representations of dimension and purpose:
typed 2D/3D containers coexist with an independent `is_3d` flag and a mutable
`source_coordinate_dim`. The goal should be fewer ambiguous states, not a
larger configuration framework.

## 2. Observed Ambiguities and Proposed Resolution

### 2.1 A `Conformer3D` can declare itself non-3D

In `crates/cosmolkit-model/src/coordinates.rs`, `Conformer3D` stores XYZ rows
and an independent `is_3d: bool`. The current validator checks row counts and
finite values, but accepts `is_3d == false`. Core stereochemical assignment
can skip such a conformer even though it is exposed through `conformers_3d()`.

Proposal: final public `Conformer3D` values should unambiguously represent
spatial conformers. Keep raw XYZ rows, original file flags, and effective
interpretation state in a private IO representation while reproducing upstream
parsing and finalization. Convert to final typed model values at a defined
boundary, without losing source-required behavior.

Do not implement this as simply deleting the flag: inventory its consumers
and preserve the relevant upstream branches first. A spatial conformer may
be planar, including having all Z values equal to zero; numerical planarity
alone does not make it a drawing layout.

### 2.2 `source_coordinate_dim` does not reliably record source provenance

`MoleculeState::try_new` in `crates/cosmolkit/src/molecule.rs` and coordinate
installation in `crates/cosmolkit/src/ops/context.rs` recompute the field from
the current collections, preferring ThreeD whenever any 3D conformer exists.
This is storage-derived state, not necessarily the dimension of the input.
One dimension value also cannot describe simultaneous 2D and 3D storage.

Proposal: decide whether the information is provenance or a derived query.

- If provenance is required, store it as import metadata with a clear scope
  and preserve it independently of later coordinate generation. Decide how
  metadata applies to multiple imported conformers.
- If only current availability matters, derive it from the collections rather
  than maintaining a redundant source-named field.

### 2.3 Multiple 2D layouts have no complete public selection contract

`CoordinateBlock` stores a vector of 2D conformers, while `coordinates_2d()`
returns the first layout's rows. This is a reasonable convenience accessor,
but it does not settle whether multiple layouts are a supported public feature.

Decisions needed:

- Whether the public model supports one primary layout or multiple selectable
  layouts.
- Whether repeated 2D generation replaces the primary layout or appends one.
- How drawing and export select a layout and preserve its metadata.

Keep `coordinates_2d()` dimension-specific and read-only. Do not add a full
layout-management API unless a concrete use case requires it. If multiple
layouts remain supported, define selection before adding consumers that each
invent their own default.

### 2.4 Conformer identity and collection position are conflated

`CoordinateBlock::remap_topology` currently renumbers conformers by collection
position while remapping atom rows. IDs are validated separately in the 2D and
3D collections. An atom-removal operation can therefore change a conformer ID
even when no conformer was removed.

Proposal: preserve conformer identity across atom-row remapping. Explicitly
distinguish ID lookup from collection indexing. Separate per-dimension ID
namespaces are acceptable when selection includes the dimension; a global ID
framework is not automatically necessary.

Document append, replacement, deletion, and ID-allocation behavior. Do not
silently reuse positional numbering as stable identity.

### 2.5 Export still makes implicit coordinate choices

At the review date, `write_v3000_detached` in
`crates/cosmolkit-io/src/sdf.rs` selected the first 3D conformer, otherwise the
first 2D layout, and otherwise synthesized zero rows. Its default dimension
header was based on collection membership rather than the conformer's flag;
a retained file-info property could also override the generated header.

Consequently, generating a 2D layout does not necessarily make a subsequent
export use that layout. This is not inherently wrong, but it needs an explicit
contract rather than an incidental implementation order.

Proposal: specify a documented default and explicit dimension/conformer
selection. Define missing-coordinate behavior and reconcile output dimension
metadata with the selected coordinates. Do not silently select another
dimension when an explicit request cannot be satisfied. Audit each writer
against its pinned source rather than treating one writer as representative
of every format.

### 2.6 Import interpretation is not output filtering or coordinate generation

The inspected historical SDF implementation mapped `Preserve` to the detected
dimension, `Require2D` to `is_3d = false`, and `Require3D` to `is_3d = true`
before molfile finalization. These flags affected stereochemical perception;
the final storage conversion happened afterwards. `Require3D` did not generate
a spatial embedding or simply reject non-3D input.

This historical observation came from matching `sdf.rs` files in sibling
`COSMolKit_1` and `COSMolKit_2` directories. Matching contents do not establish
an approved frozen baseline. Confirm provenance and authorized reference use
before relying on those snapshots for implementation acceptance.

Proposal: separate the concepts in documentation and contracts:

- Input interpretation: how file flags and coordinates participate in parsing
  and stereochemical perception.
- Layout or conformer generation: algorithms producing new coordinates for an
  already interpreted molecular structure.
- Export selection or explicit projection: which stored coordinates are
  written, or which lossy transformation was requested.

An advanced interpretation override can be useful for malformed dimension
metadata, but its name should disclose that it changes interpretation. Names
such as `Force2D` are candidates, not approved API names. Avoid implying input
validation, filtering, or 3D generation through a single ambiguous mode.

Do not replace existing semantics with post-parse filtering under the same
name. Conversely, do not preserve a confusing name merely because it existed
historically. Approve the public semantic decision and register it before
implementation.

## 3. Proposed Minimal User Contract

The following are proposals, not changes to current normative requirements.

| User intent | Proposed contract |
|---|---|
| Read SDF | Interpret input and stereochemistry according to documented source-backed rules; retain the resulting coordinates. |
| Override input dimension interpretation | Explicit advanced option whose effect on stereochemical perception is documented. |
| Generate a 2D layout | Change the selected 2D layout; preserve 3D conformers and established molecular stereochemistry. |
| Generate 3D conformers | Explicit append/replace behavior for 3D conformers; preserve 2D layouts. |
| Query coordinates | Read-only, no generation and no fallback to another dimension. |
| Export | Document default selection and allow explicit dimension/conformer selection where applicable. |
| Change topology | Apply the authoritative atom mapping to every retained coordinate set; preserve conformer IDs unless identity itself changes. |

Preserving molecular stereochemistry does not prohibit generating drawing
annotations such as wedge directions. Distinguish presentation annotations
from changing the molecule's established stereochemical identity. Explicit
stereochemical reassignment operations remain a separate capability.

Do not add independent interpretation, storage, conversion, and selection
policy objects preemptively. First settle defaults and supported use cases;
expose only the options needed to express meaningful distinctions.

## 4. Deferred Review and Acceptance Topics

When this proposal is explicitly activated, re-inspect the then-current code;
the observations above are a dated snapshot. Resolve the semantic decisions
before changing types, public signatures, registries, or ported algorithms.

Review and test at least:

- Coexisting 2D and multiple 3D coordinate sets; generation in either dimension
  preserves the other collection and its metadata.
- Planar spatial conformers, genuinely 2D depictions, inconsistent file flags,
  and explicit interpretation overrides without heuristic substitution.
- Import-time stereo perception versus later depiction generation, including
  cases where a dimension override intentionally changes source perception.
- Repeated layout generation and explicit selection, with no ambiguous default
  between stored layouts.
- Stable conformer IDs across atom removal/reordering, correct coordinate-row
  mapping, block sharing, and failure atomicity.
- Source provenance after importing and subsequently generating coordinates.
- Export selection, absent coordinates, selected-ID errors, and consistent
  dimension headers; round-trip tests must state expected information loss.
- Binding-contract registration, operation effects/invariants, exact defaults,
  and source anchors for changed behavior. Run the owning and public strict
  validation gates required by the active project rules.

This draft authorizes no immediate code changes, plan changes, commits,
releases, or instructions to the currently executing agent.
