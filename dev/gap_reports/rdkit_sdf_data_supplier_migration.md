# RDKit SDF Data And Supplier Migration

## Outcome

`cosmolkit-io` now owns the detached SDF data-field, property-list, forward
record-framing, and random-access indexing boundary:

- `SdfDataReadParams` controls strict data-field parsing and atom/bond property
  list application without claiming live-molecule chemistry policy;
- `read_sdf_*_with_params()` preserves raw molecule properties and data fields
  whether or not property lists are applied;
- `SdfGraphReader` consumes ordered concrete or query-bearing records and
  remains positioned at the next record after a parse error;
- `SdfGraphDataset` indexes byte and line spans, exposes record metadata, and
  seeks directly to an indexed record; and
- the runtime `SdfDataset` and forward `SdfReader` production paths delegate
  record indexing/framing to `cosmolkit-io` before applying runtime-owned
  MolBlock compatibility and chemistry finalization.

The detached parameter type intentionally excludes `sanitize`, `remove_hs`,
`expand_attachment_points`, and coordinate coercion. Those options require a
live `Molecule`, derived-state installation, or registered operations and
remain runtime adapter responsibilities.

## Source-Backed Behavior

The detached owner follows the pinned RDKit `ForwardSDMolSupplier`,
`SDMolSupplier`, and `FileParserUtils` paths for:

- data labels, multiline values, carriage-return removal, invalid-header
  skipping, and strict/non-strict spurious-content behavior;
- `atom.prop`, `atom.iprop`, `atom.dprop`, `atom.bprop` and the corresponding
  bond prefixes;
- default and per-list missing-value markers, item-count validation, lexical
  filtering, and boolean `1`/`0` conversion;
- delimiter recognition only when `$$$$` begins a line;
- forward recovery after a malformed record; and
- cached byte offsets, line offsets, lengths, titles, and indexed out-of-range
  errors.

Property payloads remain strings because that is the currently modeled
atom/bond property representation. Typed list intent and optional values are
also retained in `MoleculeProperties::sdf_property_lists()`.

## Performance Status

Random access seeks directly to the cached byte offset and reads one record,
matching RDKit's post-index lookup complexity. Index construction and forward
parsing remain marked performance-axis `❌`: the Rust owner reads line by line
and buffers a complete record, while RDKit's indexed supplier scans fixed-size
chunks and its classic forward supplier parses from the stream in place.

## Residual Boundaries

- `SdfDataReadParams::strict_parsing` now reaches the embedded detached
  MolBlock reader. Remaining MolBlock strict/legacy gaps are tracked in
  `rdkit_molblock_strict_legacy_migration.md`.
- The runtime still owns the compatibility MolBlock parser and the
  `finishMolProcessing` sequence. It delegates supplier framing/indexing but
  cannot delegate all parsing until the detached reader covers every legacy
  MolBlock state accepted by that public API.
- Runtime-only test helpers retain the former framing/index algorithms as
  comparison anchors; they are excluded from production builds and can be
  removed with the final runtime IO cleanup.
- Warning-log text is not modeled; invalid list values are ignored with the
  same graph-state result and remain behavior-axis `❗` where logging differs.

## Validation

```text
cargo check -p cosmolkit-io
cargo test -p cosmolkit-io --release
cargo test -p cosmolkit --release --features op-contracts-strict sdf_dataset
73 detached IO tests passed; 2 focused runtime dataset tests passed
```

Focused tests cover disabled property-list application, non-strict data-field
recovery, query-preserving forward reads, recovery after a malformed record,
metadata offsets, direct indexed reads, and out-of-range access.
