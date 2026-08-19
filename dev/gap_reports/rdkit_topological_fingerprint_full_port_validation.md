# RDKFingerprint Full-Port Validation

Status: complete for the selected RDKit `RDKFingerprintMol` bit-vector and
`atomBits`/`bitInfo` boundary. Avalon has a separate completed validation
report and is not covered by this family report.

## Reference and corpus

- RDKit: `2026.03.1`, revision `351f8f378f8ad6bbd517980c38896e66bf907af8`.
- Reference runtime: CPython `3.13.12` in `.venv/bin/python`.
- Corpus: `testdata/smiles/corpus/smiles_5000.smi`, 5,000 rows.
- Corpus SHA-256: `a4d579cd72621af27772256bb23ba796452276bb924fd20aac83625ffa67d849`.
- Golden: `testdata/fingerprint/expected/rdkit/smiles_5000/rdkit_topological_fingerprint.jsonl`.
- Golden SHA-256: `0d443196f862ad5d2fa0c9be3c16b8fc3464800d0ea64f87f894b2eba1f4cc69`.
- Golden manifest SHA-256: `83b1865dd755871afa5942df997c1d697c23e3b7a284b500625b9d6945b40651`.
- Profile: `tools/testdata/rdkit/rdkit_topological_fingerprint_profile.json`,
  14 source-meaningful parameter branches.

## Exact parity result

The generated JSONL contains 5,000 records, all with `rdkit_ok=true`, and 14
successful branches per record. The Rust test consumes every record and every
branch. It compares the source-sized parameter object, output bit-vector size,
on-bit count, and the complete sorted on-bit set. There is no sampling,
tolerance, fallback, filtered row, or aggregate similarity assertion.

Commands:

```text
RAYON_NUM_THREADS=16 COSMOLKIT_PARITY_PROFILE=smiles_5000 \
  cargo test -p cosmolkit-core --features op-contracts-strict \
  --test rdkit_topological_fingerprint_golden --release -- \
  --nocapture --test-threads=1
```

Result: 2 tests passed, including the row/branch structural check and
70,000 exact profile-row comparisons. The final run completed in 6.03s.

The maintained gate is supplemented by a full ChEMBL 37 audit against the
same pinned RDKit version. It processed 2,897,819 source records, of which
2,897,804 were mutually parseable and 15 were rejected by both libraries. All
mutually parseable records matched across the 14 vector profiles (40,569,256
comparisons) and two complete provenance profiles (5,795,608 comparisons).
The provenance profiles compare output width, complete on-bit sets,
`atomBits`, and `bitInfo`. All 128 shards completed with zero mismatch or
failed task, for 46,364,864 exact topological comparisons.
The audit ran on 2026-08-19 with COSMolKit `0.2.12` at commit `ea7c581c`;
the ChEMBL 37 source SHA-256 is
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.

Focused branch coverage also passes:

```text
cargo test -p cosmolkit-core --features op-contracts-strict \
  --test rdkit_topological_fingerprint_parity -- --nocapture
```

Result: 3 tests passed, covering exact focused bit matrices, provenance before
folding, density folding provenance preservation, and source precondition
errors. The Python API target passes 20 tests.

## Error and provenance parity

- `minPath==0`, `maxPath<minPath`, `fpSize==0`, and `nBitsPerHash==0` are
  rejected as structured invalid arguments in the focused Rust and Python
  tests.
- Query-containing paths are rejected at the source hash-input boundary and
  are covered by the focused Rust unit test; no query fallback is used.
- `atomBits` and `bitInfo` are allocated only when requested, preserve source
  atom/path ordering, deduplicate atom-to-bit entries exactly as RDKit does,
  and retain pre-fold bit identifiers.
- The 5,000-row profile has no RDKit error rows and no unexplained COSMolKit
  errors. Error text is intentionally not used as a substitute for structured
  error-kind checks.

## Source-line closure

The active source ranges are embedded in
`crates/cosmolkit-core/src/properties/fingerprint.rs` with individual
two-axis markers:

- topological helper range (`FingerprintUtil.cpp`, `Subgraphs.cpp`,
  `RDKitFPGenerator.cpp`, and Boost hash/RNG logic): 218 `RDKit✔️✔️` lines;
- public boundary and output assembly (`Fingerprints.cpp`): 25
  `RDKit✔️✔️` lines;
- `ROMol::getBondBetweenAtoms` endpoint lookup is copied at the Rust helper
  boundary and uses COSMolKit adjacency lookup, preserving source-local
  degree complexity rather than scanning every bond.

No `RDKit❌` or `RDKit❗` marker remains in the selected topological source
range. The inventory ledger in
`dev/gap_reports/rdkit_topological_fingerprint_source_inventory.md` is updated
to `✔️✔️` for every selected source unit.

## Determinism and performance review

The exact test was run twice with `RAYON_NUM_THREADS=16`; both runs passed with
6.30s and 6.60s wall time. The generated golden and manifest checksums are
stable. A fixed one-thread Rust run also passed in 64.93s. The 14-branch
single-process RDKit oracle timing was 98.61s for the same 5,000 rows.

These timings are not presented as a raw apples-to-apples speed claim: the
Rust parity target uses 16 deterministic Rayon workers while the Python oracle
is single-process. They do establish that the port has no observed
complexity regression on the maintained corpus. The source-aligned hot paths
use indexed topology access, dense atom adjacency for linear paths, local
degree bond lookup, path-local vectors, and the same recursive path-copy shape
as the pinned C++ implementation. The final output uses ordered bit storage,
matching the source's explicit bit-vector semantics.

Machine artifacts are retained under ignored `tmp/parity-audit/`:

- `rdkit-topological-oracle-smiles_5000.log`
- `rdkit-topological-oracle-smiles_5000.checksums`
- `rdkit-topological-rust-run-1.log` and `.time`
- `rdkit-topological-rust-run-2.log` and `.time`
- `rdkit-topological-rust-scalar.log` and `.time`
- `rdkit-topological-rdkit-timing.log`
- `rdkit-topological-focused-after-bond-lookup.log`

## Boundary

This report validates only the selected `RDKFingerprintMol` behavior. It does
not claim parity for LayeredFingerprint, PatternFingerprint, unfolded/count
RDK fingerprints, the public fingerprint-generator hierarchy, Avalon, or any
other unsupported fingerprint family. Avalon's separately selected explicit-
bit boundary is covered by its own completed validation report.
