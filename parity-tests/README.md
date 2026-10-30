# Parity tests

Prepare native references once, then compare COSMolKit through Cargo. Run from
the repository root with the development `.venv` installed. Reference versions
are RDKit **2026.09.1** and, for BIO, Gemmi **0.7.5**.

```bash
cargo run -p cosmolkit-parity-tests --profile dev-test --features cosmolkit/op-contracts-strict -- prepare --corpus smiles_small --threads 32
cargo test -p cosmolkit-parity-tests --profile dev-test --features cosmolkit/op-contracts-strict --test corpus -- --corpus smiles_small --threads 32
```

Both commands default to `full`. `dev-test` is optimized; `--release` also works.
Comparison never generates references. Repeating `prepare` reuses matching
references without starting Python. Use the same selection and parameters in
both commands; changing the worker count does not invalidate references.

## Select work

Append these options to **both** commands:

| Selection | Arguments |
|---|---|
| One task | `--task fingerprint_morgan` |
| Several tasks | `--task cip,valence` |
| A family | `--suite fingerprints` |
| InChI on 5,000 bundled inputs | `--suite inchi --corpus smiles_5000` |
| Several families | `--suite cip,valence,rings` |
| Fixed special regressions | `--suite special` |
| First 5,000 records of a file | `--corpus /data/chembl37.smi --limit 5000` |
| Next 5,000 records | `--corpus /data/chembl37.smi --offset 5000 --limit 5000` |
| Additional deterministic parameter samples | `--random-profiles 4 --seed 20260306` |
| Explicit complete profiles | `--profiles /path/to/profiles.json` |

List exact task and family names:

```bash
cargo run -p cosmolkit-parity-tests --profile dev-test -- list
```

A positional task name also works: `cargo test ... --test corpus -- cip`.
This target has a custom harness: options after `--` belong to the runner, not
to Cargo/libtest. Unknown names, empty selections and zero workers fail.

Profile JSON maps task names to nonempty arrays of complete parameter objects,
for example `{"valence":[{"recalculate":true}]}`. Default matrices remain
explicit in the corresponding Rust task module. Random sampling adds profiles
to those matrices; fixed special-regression recipes are not randomized.

## Inputs and results

InChI generation rejections compare the generation stage and exact engine
return code (2: error, 3: fatal), only with an empty identifier. They count as
matching source errors, not successful identifiers. Unknown failures remain
failures; successful identifiers, keys and readback state stay bit-exact.
The InChI-only `inchi-generation-status-v1` reference contract requires its
three references to be prepared once; other tasks retain their cache keys.

`testdata/corpora/*.smi` contains SMILES inputs. A bare corpus name selects that
directory; an explicit path accepts any SMILES file. Blank/comment lines are
ignored; the first whitespace-delimited field is the SMILES. Offsets and limits
count input records, not physical lines.

BIO and the special groups use their complete bundled JSONL construction
recipes under `testdata/bio/` and `testdata/special/`. They run independently of
the selected SMILES range. For a custom fixed-input selection, use exactly one
`--task` with `--corpus /path/to/cases.jsonl`; offset/limit then apply to that file.
Molecular inputs, coordinates and options belong in recipes, never expected
chemical outputs. No repository-root `testdata/`, historical handoffs or old
test-runner modules are required.

`--task tetrahedral_geometry` restores the historical ordered-ligand
geometry check. Prepare uses RDKit ETKDGv3 with seed 42; comparison checks CK
center indices and `det(ligand0-center, ligand1-center, ligand2-center) > 1e-8`
against those independent coordinates, through scalar and batch APIs. Its 30
inputs exclude `historical_23`, whose reference embedding failed. This is
a handedness predicate, not a tolerance on coordinate equality. Failed reference
embedding remains a reported unavailable comparison, never a pass.

Generated references live in ignored `expected/v1/<identity>/`; streamed reports
live in ignored `reports/run-*/`. Summaries distinguish matching values,
matching source-defined errors, differences, native crashes, undefined outcomes
and timeouts. Worker timeouts and embedding's explicit `EXCEEDED_TIMEOUT` counts
on either side are recorded separately: neither failures nor matches, even when
both sides time out. They do not cause a failing exit; a timeout-only task is
marked `no_comparable_results`, not a parity pass. Other exceptional native
outcomes are also not parity passes. Reports retain the input, full profile and
first differing byte; large payloads retain hashes and bounded byte windows.
The reference file retains the complete expected result for reproduction.

Every numeric field is compared by its exact encoded IEEE bits, including
coordinates, energies, gradients and alignment matrices. There is no epsilon,
rounding or numeric hash-only acceptance. Exception wrapping may differ only
for explicitly paired native diagnostics and concrete Rust error kinds.
Declared representation projections (such as MACCS's unused zero bit and
typed property slots) are explicit in adapters/profiles, not chemical fitting.

The `inchi` suite compares exact InChI strings and both molecule-to-key and
identifier-to-key calls with default, `/FixedH` and `/SNon` options. Its
roundtrip task compares the generated identifier and every atom/bond field
used by the core state adapter after reading it with all four
`sanitize`/`remove_hs` combinations. This is a SMILES → InChI → molecule chain,
not an independent parser test using a shared external InChI input. Native
InChI calls use the existing crash/timeout isolation.
RDKit's Python reader turns sanitization exceptions into `None`; this recorded
read rejection matches CK's typed `mol_from_inchi` `SanitizeFailed` outcome.
Both diagnostics remain in the report. Internal exception fields are not
compared when the reference API does not expose them; unrelated errors do not match.

MMFF optimizers run directly on the current state. Native MMFF property
construction mutates its molecule, so energy observations prepare a temporary
copy and cannot precondition the next optimization. The optimizer profiles
record this as `optimizer_observation: isolated_queries_v1`.
UFF single optimization obtains status and energy from the same native force
field, as `FFConvenience.h::OptimizeMolecule` does; rebuilding it after minimize
is a different observation. Its profile records
`optimizer_observation: retained_optimizer_v1`.

## Execution and cache

Tasks run **serially**. `--threads N` selects N persistent native **processes**
for preparation, avoiding Python's GIL, and N Rayon workers for comparison.
Equivalent public batch APIs receive the same N: the first ten corpus inputs
use scalar APIs, and the remainder use batch. Tasks without equivalent batch
APIs run inside the selected Rayon pool. Numerical calls avoid multiplying the
outer width with uncontrolled inner pools. The 77 environment-sensitive stereo
recipes isolate each Rust call in the already-built executable; no nested Cargo.

The cache key contains the pinned reference version, task, complete parameter
profiles and a checksum/count of the actual selected inputs. BIO profiles also
include Gemmi's version. It does not hash unrelated source files or depend on
thread count, timestamps or environment variables. An incompatible adapter
output change requires a new explicit result/profile identity.

Warm preflight hashes inputs once per distinct selection and reads only small
manifests and reference headers. Compressed frame checksums and the final file
checksum are verified during the one streaming comparison pass, not by loading
or rehashing the full reference first. Outputs use compact typed binary records
and fast zlib compression. Frame/record allocation limits are 128 MiB; windows
and pending jobs are bounded independently of corpus length. Native workers
compress to temporary files instead of accumulating expanded results in IPC
queues. Inputs/profiles are not expanded into millions of JSON objects.

## Framework regressions

```bash
cargo test -p cosmolkit-parity-tests --profile dev-test --features cosmolkit/op-contracts-strict --lib
.venv/bin/python -m unittest discover -s parity-tests/tools -p 'test_*.py'
```

These check selection, reproducibility, exact encoding, streaming checksums,
allocation limits, fixed-case censuses and native adapter entrypoints. Ordinary
chemistry regressions remain in their owning crates.
