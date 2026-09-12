# RDKit SDF Typed SGroup Migration

## Outcome

`cosmolkit-io` now owns detached V2000 and V3000 substance-group parsing and
writing over `cosmolkit-model::SubstanceGroup`. V3000 enhanced-stereo
collections are parsed into and written from `TopologyBlock::stereo_groups`.
These records no longer fail at the concrete detached boundary or depend on
the legacy runtime SDF parser merely to retain typed topology annotations.

The migrated reader covers:

- V3000 `ATOMS`, `PATOMS`, `XBONDS`, `CBONDS`, `BRKXYZ`, `CSTATE`, `SAP`,
  `PARENT`, `COMPNO`, data fields, string properties, defaults, non-sequential
  atom/bond bookmarks, and unsorted SGroup sequence IDs;
- V3000 `STEABS`, `STEREL`, and `STERAC` collection records, including the
  strict single-absolute-group rule;
- V2000 `STY`, `SST`, `SLB`, `SCN`, `SDS`, `SAL`, `SBL`, `SPA`, `SMT`, `SDI`,
  `SBV`, `SDT`, `SDD`, `SCD`, `SED`, `SPL`, `SNC`, `SAP`, `SCL`, and `SBT`
  property records; and
- delayed parent resolution after sequence-ID sorting.

The detached writer emits the corresponding V2000 property records and V3000
SGroup/collection blocks. It preserves contained versus crossing V3000 bond
roles, assigns missing or duplicate enhanced-stereo write IDs by group kind,
uses RDKit's continuation threshold, and wraps V2000 UTF-8 data fields without
splitting a scalar value.

## Source Status

The new implementation is aligned against pinned RDKit
`MolSGroupParsing.cpp`, `MolSGroupWriting.cpp`, and the enhanced-stereo paths in
`MolFileParser.cpp` and `MolFileWriter.cpp`. Corresponding Rust functions carry
verbatim source excerpts with two-axis markers. Most aggregate SGroup paths
remain `RDKit❗✔️`: the typed state and tested strict behavior are present, but
non-strict warning-and-drop recovery, every arbitrary extension property, and
atropisomer-derived collection atoms have not been proven equivalent.

Unsupported query-bearing atom/bond records and V3000 `OBJ3D` constraints
remain structured errors. They are separate from SGroup support and must not
be approximated into concrete topology.

## Validation

```text
cargo check -p cosmolkit-io
passed

cargo test -p cosmolkit-io --release
49 passed, 0 failed
```

Focused tests cover V2000 and V3000 typed read/write/read behavior, hierarchy,
non-sequential bookmarks, unsorted sequence IDs, enhanced stereo, declared
SGroup counts, and a multibyte UTF-8 V2000 data-field wrapping boundary.

## Remaining Ownership Work

The legacy runtime MolBlock/SDF parser still owns parameterized sanitization,
hydrogen removal, streaming/indexed suppliers, property-list processing, and
the broader source parser surface. It must remain until the detached owner has
an explicit query-bearing MolBlock result contract and the facade can apply
finalization policy without retaining a second parser.
