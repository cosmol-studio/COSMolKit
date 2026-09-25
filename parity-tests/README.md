# Rust parity registry pilot

This top-level, unpublished crate tests the **public `cosmolkit` API with
`full` enabled**. Its only chemistry dependency is `cosmolkit`. No production
crate depends on this runner. Scope is exactly `fuzzy_and` and `fuzzy_or` on
the two public sparse-count index widths; this is not a performance suite.

## Read the implementation in this order

1. `src/registry.rs`: two tasks, their u32/u64 widths, typed inputs/results,
   and six explicitly named branch cases. No YAML and no hidden repetitions.
2. `src/lib.rs`: `run` automatically prepares data, performs complete preflight,
   then executes and compares exact typed results.
3. `src/execute.rs`: calls the public Rust API; no chemistry reimplementation.
4. `src/tests.rs`: tests the runner itself, including missing final-task data
   preventing **all** Rust operation calls.

The RDKit adapter in `tools/oracles/rdkit/fingerprint_values_pilot.py` calls
RDKit 2026.03.1 only. It is not a COSMolKit Python binding test. Rust selects
and expands tasks; Python receives the exact typed rows and invokes RDKit.
Python/JS binding comparison is deliberately not implemented in this pilot.

## Run from the repository root

```bash
# Inspect the entire matrix; no oracle or CK operation runs.
cargo run -p cosmolkit-parity-tests --release -- list

# Optional cache preparation without CK execution. Not a prerequisite for run.
cargo run -p cosmolkit-parity-tests --release -- prepare

# Validate all selected artifacts without executing CK operations.
cargo run -p cosmolkit-parity-tests --release -- preflight

# Prepare/reuse ALL tasks, globally preflight, then execute Rust and compare.
cargo run -p cosmolkit-parity-tests --release -- run

# Focus on one task, for preparation or execution.
cargo run -p cosmolkit-parity-tests --release -- prepare --task fuzzy_and
cargo run -p cosmolkit-parity-tests --release -- run --task fuzzy_and

# Runner regressions use synthetic references; no RDKit/environment required.
cargo test -p cosmolkit-parity-tests --release
```

Use `--python PATH` on run or prepare to select a reference interpreter. Its RDKit
version must match the registry. Generated inputs, reference results,
manifests and Rust results default to ignored `target/parity-tests/`.
`--data DIR` selects an external/ignored artifact directory, not a Git input
directory. Do not stage these generated artifacts.

## What exactly executes?

The builtin input pairs exercise:

| Case | Source branch |
|---|---|
| `empty_both` | Neither loop has entries |
| `left_empty_right_tail` | Union inserts right tail; intersection stays empty |
| `right_empty_remove_left` | Intersection removes left; union preserves it |
| `shared_signed_min_max` | Shared keys compare signed values, with exclusive tails |
| `disjoint_interleaved` | Ordered scan advances across alternating keys |
| `explicit_zero_shared_and_exclusive` | Stored zero is an entry, not an absent key |

Each case executes once for each registered width and operation:
6 input pairs x 2 widths x 2 operations = **24 comparisons**. Selecting one
operation gives 12. These counts follow the declared cases, not a target
sample count, repetition constant, or performance claim.

Supply `--corpus FILE` to all stages to use another complete input list.
The file is a JSON array of the Rust `Pair` type; operation/width expansion
still comes from the Rust registry, not the file:

```json
[
  {
    "id": "signed_overlap",
    "length": 16,
    "left": [[1, 5], [3, -2], [8, 4]],
    "right": [[1, 3], [3, -4], [9, 7]]
  }
]
```

This pilot accepts equal-length valid vectors with unique indices strictly
below length and i32 counts. Construction creates stored zeros separately
before assigning nonzero counts; no arbitrary count-magnitude limit is used.
Other inputs fail preflight rather than being skipped. Exhaustive integer-domain coverage, error parity,
SMILES-to-fingerprint conversion and million-row streaming are not claimed.
The small pilot loads owned snapshots in memory to keep the execution flow
inspectable; it is not yet the large-corpus engine.

## Fail before doing expensive work

Preparation writes each generation into a new directory and publishes it by
rename. It records SHA-256 identities for input, reference output, Rust
registry and oracle adapter, plus version, schema and row count. A partial
generation cannot appear ready. `run` and `prepare` first inspect all selected
tasks, reuse verified generations without invoking Python, and generate missing
or invalid ones using the pinned oracle. Invalid generations are moved to
`.invalid-*` directories with their paths reported, preserving failure evidence.
A filesystem lock serializes publishers. Temporary directories prevent partial
writes from appearing ready. A comparison failure never regenerates references.

Preflight checks **every selected task** before returning executable work:

- all input/reference/manifest files exist;
- identities, exact input bytes, row counts and case/operation/width match;
- results contain valid lengths and sorted unique indices;
- no empty corpus, duplicate case IDs or invalid values are accepted.

`preflight` alone remains read-only and rejects missing data. `run` prepares it
automatically; failed generation or global preflight prevents **all** CK calls.
No task starts comparison while another task still needs preparation. Inputs/reference
records are owned snapshots after preflight. Runtime comparison collects all
selected results, including errors and mismatches; none are skipped. Output
is `rust-report.json` with typed input, expected, actual and match status.
No runtime result cache is reused. Task/variant identity remains attached to
each result so a future binding adapter can consume precisely the same cases.

## Separate fixed regressions

Fuzzy testing has exactly two categories:

1. Fixed regressions: literal expected values, errors, stored-zero and signed
   count behavior, unchanged operands, and a fixed cross-node retain case.
   These run in the owning crate; public API wiring and runner mechanics have
   their own small fixed regressions at their respective boundaries.
2. Registry-driven parity: corpus expansion, reference preparation and RDKit
   comparison belong exclusively to this top-level crate and its oracle adapter.

There is no separate Fuzzy performance sampler or owner-local size/order/pattern
corpus sweep. The former 540-case sweep was replaced by one fixed retain
regression per width. Tests formerly named `performance_fuzzy_*_shapes_both_widths`
are named `regression_fuzzy_*_fixed_shapes_both_widths`: they assert fixed results,
not performance. Historical audit logs do not define additional test entrypoints.

`crates/cosmolkit/tests/fingerprint_values_api.rs` checks public export and
contract wiring with fixed input/expected values. Owner regression tests stay
in `crates/cosmolkit-fingerprints/tests/`. Neither contains this pilot's corpus
loader, preparation, oracle invocation or comparison pipeline. The old
1400-call timing sampler was removed; its functional shapes remain covered
by deterministic owner regressions. No speed equivalence is inferred.

## Public API audit

`binding_contract/registry.rs` registers both types and every exposed method
before exposure. These are direct re-exports of the owner types, not wrapper
values; private index validation remains private. No old-name aliases exist.

| Canonical Rust method | Pinned RDKit public source |
|---|---|
| `new` | `SparseIntVect(length)` |
| `length`, `value`, `set_value`, `nonzero_elements`, `total_value` | `getLength`, `getVal`, `setVal`, `getNonzeroElements`, `getTotalVal` |
| `fuzzy_and`, `fuzzy_or` | non-mutating `operator&`, `operator\|` |
| `with_added`, `with_subtracted` | vector `operator+`, `operator-` |
| `with_added_scalar`, `with_subtracted_scalar`, `with_multiplied_scalar`, `with_divided_scalar` | value-returning CK projections of public scalar `+=`, `-=`, `*=`, `/=` |

Source files: `Code/DataStructs/SparseIntVect.h` and
`Code/DataStructs/Wrap/SparseIntVect.cpp:135-185` in the pinned RDKit tree.
The existing approved structured errors for source undefined arithmetic are
unchanged. Detached value setters do not acquire Molecule commit authority.
The registry's type-owned `read_only` state describes the borrowed receiver,
including methods returning a new value; no macro checker was weakened.
Python/JS names in the binding contract are declarations, not proof that
those language adapters are implemented.

## Automatic preparation validation

- Release runner regressions: 13 passed, 0 failed/ignored (exit 0).
- Cold `run` with RDKit: 2 tasks generated, 24 comparisons, 0 failures (exit 0).
- Warm `run` with a nonexistent Python path: 2 tasks reused, no oracle calls,
  24 comparisons, 0 failures (exit 0).
- Cold `run` with a nonexistent Python path: preparation failed, 0 CK calls
  (expected exit 1).
- Regression coverage includes invalid-reference preservation and repair,
  stale manifests, malformed oracle output, selected-task isolation, the
  all-task execution barrier, and no reference rewrite after CK failures.

## Earlier API validation (before automatic preparation)

- `cargo test -p cosmolkit --release --features full,op-contracts-strict`:
  297 passed across the executed test targets (including 3 compile-fail
  doctests); 0 failed/ignored. Zero-case targets are not counted as coverage.
- Owner integration regressions: 104 passed, 0 failed/ignored after removing
  the timing sampler. The standalone target is
  `cargo test -p cosmolkit-fingerprints --release --test migration_fp_values`.
- `cargo check -p cosmolkit --features full,op-contracts-strict`: exit 0.
- `cargo check -p cosmolkit --no-default-features`: exit 0.
- Default `prepare`, `preflight`, `run`: 24 RDKit comparisons, 0 mismatches.
- Standalone preflight rejects a missing fuzzy_or generation before CK calls.
  The automatic `run` pipeline now prepares that missing generation instead.
- `cargo fmt --all -- --check`: exit 0.

The earlier exact registry-order test failed because its expected list lacked
the 31 new entries. The expected list was extended explicitly without removing
its exact order/count assertions; the focused 7-test schema target and full
runtime rerun then passed. Existing compiler warnings remain. These results
are not a whole-workspace, Python/JS, performance or million-corpus claim.
