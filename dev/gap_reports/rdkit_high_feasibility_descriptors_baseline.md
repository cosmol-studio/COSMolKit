# High-Feasibility Descriptor Port Baseline

## Command

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict
```

The baseline was captured before any descriptor implementation or core
refactoring from this plan.

## Result

```text
3515 passed
5 failed
46 ignored
0 measured
0 filtered out
```

The release build completed in 5 minutes 37 seconds. The library test binary
completed in 1.28 seconds and exited with status 101 because of the five
failures below.

## Pre-Existing Failures

| Test | Owner | Failure |
|---|---|---|
| `properties::descriptors::qed_tests::delete_substructs_golden_covers_required_parameter_branches` | Substructure expected data | Unified generator checksum changed from the checksum recorded in the committed expected-data manifest. |
| `properties::descriptors::qed_tests::delete_substructs_matches_rdkit_golden_for_only_frags_matrix` | Substructure expected data | Same generator-checksum mismatch. |
| `properties::rdkit_prepared_draw_parity::prepared_draw_golden_has_one_record_per_smiles` | Depiction expected data | Same generator-checksum mismatch. |
| `properties::rdkit_prepared_draw_parity::prepared_draw_molecule_matches_rdkit_golden` | Depiction expected data | Same generator-checksum mismatch. |
| `properties::rdkit_prepared_draw_parity::prepared_draw_molecule_matches_rdkit_golden_in_parallel_batch` | Depiction expected data | Same generator-checksum mismatch. |

Observed generator checksum:

```text
2469e2afc6694eb2579b87679d67dec4704f20b6d2e4b8064d6b562e0a39ab4c
```

Checksum recorded by the affected expected data:

```text
8e2bb795f66626be98d89b8e592a382689bed107055f1bf59d064e044c5e25e4
```

## Disposition

These failures predate the descriptor port and do not indicate a descriptor,
substructure, or depiction output mismatch. They are reproducibility-gate
failures caused by stale expected-data provenance after the unified generator
changed on the remote baseline.

The descriptor plan must not weaken or bypass this gate. The owning expected
data should be regenerated through the commands reported by the tests:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py --python .venv/bin/python --profile smiles_small --suite substructure --jobs 4
.venv/bin/python tools/testdata/rdkit/generate_all.py --python .venv/bin/python --profile smiles_small --suite depiction --jobs 4
```

Final descriptor signoff requires either a clean strict baseline after that
regeneration or an explicitly retained external blocker with proof that the
descriptor changes introduced no additional failure.

The build also emitted existing compiler warnings across the source-reproduced
InChI and core code. They are not test failures and are outside this plan, but
the descriptor port must not add new warnings.
