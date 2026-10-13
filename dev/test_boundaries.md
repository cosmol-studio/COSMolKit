# Regression Tests and Corpus Parity — 0.5.0

There are two test categories:

- Ordinary regressions verify local behavior, boundaries and errors in the
  owning crate's unit modules or `tests/`. Python and JavaScript tests focus
  on their language boundary rather than repeating complete chemistry corpora.
- Corpus tests and explicitly designated special regressions use the current
  `parity-tests/` workflow.

Completed validation and remaining coverage are recorded in
[VALIDATION.md](../VALIDATION.md).

## Ordinary regressions

Use small fixed inputs and source-backed expected values directly in tests.
Larger shared fixed fixtures may remain under repository `testdata/` with
provenance, licenses and checksums. Pinned `third_party/` fixtures may be read
at runtime under the repository policy, but production/package builds must
not require them. Do not embed external checkout paths at compile time.

Ordinary tests do not generate references, invoke or compile reference
implementations, parse upstream source into expectations, or require prepared
corpus caches. Missing required fixtures are errors, not skips. A fixed table
is not a corpus merely because it has many rows.

Check complete results and required preservation: exact atom/bond state,
coordinate bits where required, typed errors, ordering, mappings, COW sharing
and failure atomicity. Successful parsing alone does not establish a roundtrip.
Keep discovered counterexamples at the smallest owning behavior boundary.

Compiler-level private-module isolation retains its real-layout compile-pass
and compile-fail checks as required by the operation standard; runtime
rejection or textual matching cannot replace that boundary. Other ordinary
regressions use Cargo's normal test scheduling without nested builds.

## Corpus tests

The current runner's executable registry defines tasks and comparisons.
Corpus inputs belong in `parity-tests/testdata/`, generated references
in `expected/`, and reports in `reports/` under the same package.

From the repository root:

```bash
cargo run -p cosmolkit-parity-tests --profile dev-test --features cosmolkit/op-contracts-strict -- prepare --corpus smiles_5000 --threads 112
cargo test -p cosmolkit-parity-tests --profile dev-test --features cosmolkit/op-contracts-strict --test corpus -- --corpus smiles_5000 --threads 112
```

Use the same `--task NAME` or `--suite NAME` selection in both commands.
Comparison options follow Cargo's `--` separator and belong to the custom
runner, not to Cargo/libtest. Both commands default to `full`.

Preparation runs pinned reference adapters, supports configurable worker
counts, and saves ordered inputs, parameters, reference outputs and identities.
Tests read validated snapshots and never generate or overwrite expectations.
Validate the input identities and every selected reference header before the first CK call;
preflight is read-only and cannot start reference generation.
Verify compressed-frame and final-file checksums during the streaming comparison.
Missing, corrupt or stale references fail explicitly, never by regenerating data
during comparison. Tasks run serially; cases within each task use the selected
Rayon pool or equivalent public batch API. Numeric outputs compare by exact
IEEE bits; do not relax comparisons or remove cases to hide errors.
Matching source-defined errors is checked under each task's explicit error
contract, not silently counted as a successful numerical calculation.

## Special regressions

Special regressions use fixed, explicitly selected reference matrices rather
than expansion over a SMILES corpus. They have the same two-stage shape:

```bash
cargo run -p cosmolkit-parity-tests --profile dev-test --features cosmolkit/op-contracts-strict -- prepare --suite special --threads 112
cargo test -p cosmolkit-parity-tests --profile dev-test --features cosmolkit/op-contracts-strict --test corpus -- --suite special --threads 112
```

The runner documents the supported selectors, including the 77-case
`structure_tags` matrix, `tautomer_long_conjugated`, `tautomer_focused`,
`molalign_focused` and `bio_mmcif_switches`. Preserve their fixed case census,
source error branches and complete comparison fields. Test-only detached
algorithm dependencies may preserve these designated source boundaries;
this does not authorize exposing runtime internals or a second chemistry API.

Preparation and comparison instructions live only in the
[runner README](../parity-tests/README.md). Domain fixture READMEs record
provenance and purpose, not separate preparation pipelines. A missing reference
never authorizes refreshing expectations from CK results. Known failures
remain executable and are not passing evidence. Registration, zero matches,
ignored tests and historical receipts do not establish current validation.
