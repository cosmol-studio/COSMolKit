# Avalon Fingerprint Full-Port Validation

This report closes the active Avalon validation gate for the pinned source
boundary in `dev/plans/rdkit_topological_avalon_fingerprint_port_plan.md`.

## Reference boundary

- RDKit: `2026.03.1`
- Public adapter: `rdkit.Avalon.pyAvalonTools.GetAvalonFP`
- Avalon engine: `AvalonToolkit_2.0.5-pre.3`
- Engine provenance and license: [`avalon_fingerprint_source_provenance.md`](avalon_fingerprint_source_provenance.md)
- Source inventory and call graph: [`avalon_fingerprint_source_inventory.md`](avalon_fingerprint_source_inventory.md)

The generated oracle intentionally uses the pinned RDKit Python adapter. A
direct `SMIToMOL` harness remains available as a low-level engine diagnostic,
but is not used as the public parity oracle because RDKit's `ROMol` to REACCS
MolBlock conversion is the selected API boundary. This distinction matters
for query, isotope, charge, and stereochemical inputs.

## Exact corpus result

Command:

```text
RUSTFLAGS=-Awarnings cargo test -p cosmolkit-core --release \
  --features op-contracts-strict \
  --test rdkit_avalon_fingerprint_golden \
  avalon_fingerprint_matches_every_active_5000_row_profile_exactly -- --exact
```

Result: `1 passed`, 5000 corpus rows, 23 active profiles, 115,000 exact
bit-vector comparisons, zero unexplained mismatches. The test executes every
row and every profile, including query mode, all individual flag families,
combined profiles, non-SSS profiles, and the source-defined rejection rows.
It uses no sampling, tolerance, fallback, similarity threshold, or mismatch
filter. Release execution completed in 171.46 seconds on the validation host;
the deterministic parallel rerun completed in 12.97 seconds after compilation
and produced the same result.

A three-run, one-profile Python API comparison over all 5,000 rows measured a
median of 4.696 seconds for pinned RDKit and 8.411 seconds for COSMolKit. Both
paths include their public Python conversion/call overhead. The Rust port is
1.79x slower on this profile, with the same exact output; no complexity-class
regression was observed, but this constant-factor gap remains visible rather
than being described as performance parity.

Focused boundary command:

```text
RUSTFLAGS=-Awarnings cargo test -p cosmolkit-core --release \
  --features op-contracts-strict \
  --test rdkit_avalon_fingerprint_golden \
  avalon_fingerprint_matches_focused_profiles_and_size_boundaries -- --exact
```

Result: `1 passed`, covering Python/C++ defaults, query mode, non-SSS-only
flags, and 9/31/33-bit byte-rounding boundaries.

## Oracle artifacts

The generator and profile are committed under `tools/testdata/rdkit/`; the
machine-sized expected data remains ignored under `tmp/parity-audit/`.

- `avalon_fingerprint_smiles_5000.jsonl`: 5,000 rows × 23 profiles,
  SHA-256 `235b66b30ebb73e966c540e4d89f74736baaa72b4ef2a25e0d592e5d77ba7c2f`
- `avalon_fingerprint_focused.json`: 8 focused cases,
  SHA-256 `12b131ed3e1722652f6d9a274f2193727fe9c001fb6f47d52e06426ad983adfb`
- Generation metadata and command: `tmp/parity-audit/avalon_fingerprint_generation.log`

The focused and corpus records contain exact `on_bits`, requested width,
query mode, flag mask, and conversion-error status. Non-byte-aligned widths
retain the requested public width while preserving the source `nBits / 8`
byte allocation behavior.

## Source and implementation audit

- Rust Avalon modules contain individual two-axis `Avalon` source markers for
  the adapter, REACCS conversion, preprocessing, traversal, hashing, rings,
  aromaticity, and every active flag family.
- The implementation preserves source ordering, integer-width behavior,
  four-byte internal rounding, query/non-query second-pass behavior, and
  non-SSS query gating.
- `Molecule::avalon_fingerprint` and the Python wrapper are value-style and
  leave the input molecule unchanged; focused Rust/Python tests cover repeated
  calls and non-mutation.
- The source archive's BSD-like Avalon notice and RDKit adapter license are
  recorded in the provenance report; no source archive is redistributed from
  the ignored oracle directory.

## Remaining scope

Avalon count-fingerprint and string-input overloads, coordinate generation,
CheckMol, and unrelated Avalon APIs remain explicitly out of scope. They are
not implied by the explicit-bit support status below.
