# ChEMBL 37 Parity Audit

This directory is the repository-owned entrypoint for COSMolKit's large
ChEMBL 37 differential parity audit. It replaces the original ad hoc scripts
under ignored `tmp/parity-audit/` without committing the externally
distributed corpus or
multi-gigabyte run output.

The audit compares the installed COSMolKit Python extension with pinned RDKit
`2026.03.1`. It is release evidence, not an ordinary per-commit test. The
committed 152- and 5,000-record suites remain the fast and maintained gates.

## Prerequisites

- CPython `3.13.12` with COSMolKit built from the checkout under audit. The
  interpreter pin is required because exact-bit QED parity includes the
  interpreter's floating-point reduction behavior.
- RDKit `2026.03.1`, NumPy, and the Avalon RDKit module in the same environment.
- The official `chembl_37_chemreps.txt.gz` file. Download it from the URL
  recorded by `prepare_corpus.py`; verify that its SHA-256 is
  `ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.
- Enough space outside tracked repository data for roughly one gigabyte of
  prepared shards plus audit logs and summaries.

The runner uses `sys.executable`. Activate the intended environment or invoke
the runner with that environment's Python. It refuses Python or RDKit versions
other than the profile pins and refuses a COSMolKit extension whose version
differs from the workspace version. The installed NumPy version is retained in
the run identity.

## Prepare The Corpus

Preparation is deterministic and atomic. It validates the pinned source,
requires the four expected TSV columns, retains every row with a nonempty
canonical SMILES, and assigns rows by
`blake2b-64(chembl_id) modulo shard_count`.

```bash
python dev/tools/chembl_parity/prepare_corpus.py \
  --input /data/chembl_37_chemreps.txt.gz \
  --output-dir /data/cosmolkit/chembl37-128 \
  --shards 128
```

The output `manifest.json` records source identity, selection, counts, and the
size and SHA-256 of every shard. The source archive, shards, and run results
must remain outside Git; `tmp/` is acceptable for disposable local runs.

## Run The Complete Profile

Build/install the current checkout first. One supported local route is:

```bash
maturin develop --release --manifest-path python/Cargo.toml \
  --features release-abi3-py39
```

Then run:

```bash
python dev/tools/chembl_parity/run.py \
  --corpus-dir /data/cosmolkit/chembl37-128 \
  --run-dir /data/cosmolkit/runs/$(git rev-parse --short HEAD) \
  --workers 96
```

Run the command from the repository root. Worker processes are launched as
`dev.tools.chembl_parity` modules so their shared comparison helpers use
normal package imports rather than depending on the current script directory.

`profiles/complete.json` is the default 34-phase profile. It retains the 22
phases from the original completed audit, the later full-corpus modern
CIPLabeler phase, the RDKFingerprint/Avalon, AtomPair, Pattern, Topological
Torsion, and Layered fingerprint phases, and the ordinary MolAlign and tautomer
enumeration phases. Four additional phases cover topology operations,
connected fragments, direct 2D coordinate generation, and COSMolKit binary
roundtrips. The complete profile therefore contains 4,352 ordered shard tasks
over the 128-shard corpus.
The existing batch phase also compares scalar and batch sanitize, both
kekulize branches, and direct 2D coordinates. The profile records phase order,
selection boundaries, branch modes, fixed repeat phases, reference version,
expected processed counts, and the only two informational mismatch keys.
Those two keys compare historical InChI fields stored by ChEMBL with current
pinned RDKit output; COSMolKit is still required to equal RDKit on the
corresponding operation.

The four additional phases received a complete configured-corpus execution in
the final 15-phase validation run. Fragments and direct 2D coordinates had
zero blocking mismatches. Topology operations exposed unsanitized `RemoveHs`
valence/H-state differences, and binary roundtrips exposed
post-deserialization hash/fingerprint differences. The immutable identity,
exact counters, and acceptance consequences of that discovery execution are
recorded in
[`dev/gap_reports/chembl37_complete_validation_summary_20260820.md`](../../gap_reports/chembl37_complete_validation_summary_20260820.md).
Any mismatch remains blocking until the implementation is source-traced,
fixed, regression-tested, and the complete affected phase is rerun.

Both gaps completed that process. The topology implementation now reproduces
the upstream non-strict property-cache update before hydrogen removal and
migrates every surviving atom's retained valence/H state through the topology
mapping. Its complete affected-phase rerun evaluated 2,854,376 records across
all 128 shards and recorded 45,669,848 matching observations with zero
mismatch across the value-style and in-place operation branches.

The archive now retains validated derived chemistry state in a required
versioned section and remains able to read direct legacy payloads and archive
v1.0. Its complete affected-phase rerun evaluated 524,288 eligible records
across all 128 shards and recorded 11,534,336 matching observations with zero
mismatch through both public deserialization entrypoints. The binary phase is
a COSMolKit roundtrip invariant over ChEMBL structures; it is deliberately not
described as RDKit binary-format parity. These complete reruns establish the
current accepted results without rewriting the immutable discovery-run
identity above.

The topology phase deliberately selects from unsanitized parses. All 2,897,819
source SMILES are accepted by both raw parsers; 2,854,376 have at most 80 atoms
and enter the phase. This includes 14 records beyond the ordinary sanitized
`<=80` boundary so raw-success/sanitize-rejection behavior is retained.

Use one or more `--phase NAME` arguments for a preflight or a focused full
corpus rerun. A selected phase still processes every configured shard and
retains its complete branch matrix.

## Resume And Validate

Interrupted work can be resumed only with an identical run identity:

```bash
python dev/tools/chembl_parity/run.py \
  --corpus-dir /data/cosmolkit/chembl37-128 \
  --run-dir /data/cosmolkit/runs/COMMIT \
  --workers 96 \
  --resume
```

The identity includes the Git commit and status, a retained tracked-diff
patch, installed extension, profile, complete audit script set, fingerprint
profiles, corpus manifest, Python, COSMolKit, and RDKit. A task is reusable
only when its command and output checksum still match. Every corpus shard is
checksummed before both a new run and a resumed run.

Each task writes a log, JSON summary, and bounded JSONL finding examples.
Each phase writes `aggregate.json`; the run writes `manifest.json`. The runner
returns nonzero when any task is missing or fails, or when any noninformational
`mismatch.*` counter is nonzero. A time limit is never treated as completion.

## Accepting Evidence

A result may update `VALIDATION.md` only when:

1. the run manifest has `complete: true` and `passed: true`;
2. all configured phases and shards completed;
3. corpus, profile, scripts, Git state, extension, and reference identities
   are retained in the manifest;
4. any discovered mismatch has been source-traced, fixed, added as a focused
   regression, and rerun over its complete affected branch boundary; and
5. the documentation distinguishes a retained-case replay from a subsequent
   full-corpus rerun.

The profile is an explicit claimed boundary, not a mechanism for filtering
failing molecules. Do not add molecule-specific exclusions, broad exception
allowances, or new informational mismatch keys to obtain a passing run.
