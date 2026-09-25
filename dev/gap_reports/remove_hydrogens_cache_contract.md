# RemoveHs final valence cache contract

Decision: CK-VALENCE-001, explicitly approved on 2026-09-22.

`sanitize=false` preserves the complete non-sanitize removal algorithm,
including intermediate valence calculations needed for stereo and removal.
It does not run a final whole-topology calculation merely to populate a cache.
The detached result carries `final_valence: None`; the registered operation
clears both the persistent VALENCE validity bit and payload before atomic commit.

For `sanitize=true`, core computes a complete assignment after all topology and
chiral-H changes. A calculation failure propagates as an error, never None.
The operation installs a supplied final assignment; if none is available it
clears the cache. Valid means usable for the current topology under the
assignment's calculation contract, not proof of complete strict sanitization.

This intentionally differs from RDKit's observable stale pre-removal cache.
Pinned RDKit 2026.03.1, explicit `[H]C` parsed with removeHs=false: carbon's
explicit-valence/implicit-H values are 1/3 before and after
RemoveHs(sanitize=false), then 0/4 after UpdatePropertyCache(false).
CK must not install the 1/3 source snapshot as a valid final-topology cache.
This decision does not authorize changing removal topology, stereo, query
semantics, error handling, or any other source behavior. It is separate from
the degree-one coordinate exception CK-COORD-001.

Implementation scope: core RemoveHydrogensResult and its final assignment;
the thin registered hydrogen operation; focused core, public and private-cache
regressions. No new runtime primitive, public cache getter or mutable state
escape is introduced. The existing operation_defined permission already allows
update or clear; its allow-list and generic machinery are left unchanged in
this bounded implementation, not expanded. Other sanitize=false constructors
and a general effect-category cleanup are not implemented by this change.

Coverage retains the previous default/sanitize, chiral-H normalization and
failure-atomicity cases. New private-runtime assertions check both validity
and payload, source cache sharing, repeat non-sanitize deletion, explicit later
recomputation and in-place/value equivalence. These assertions do not substitute
for real-module operation privacy tests.

## Validation

Actual process results (2026-09-22; no inferred exit codes):

- Core focused: `cargo test -p cosmolkit-core --release --features op-contracts-strict --test migration_h_remove_state`, exit 0, 12 passed/0 failed/0 ignored; `/tmp/ck-valence-core-target.log`.
- Core full: `cargo test -p cosmolkit-core --release --features op-contracts-strict`, exit 0, 466 passed/0 failed/0 ignored across 35 target summaries; zero-test targets are not behavioral coverage. Log `/tmp/ck-valence-core-full.log`.
- Core check: `cargo check -p cosmolkit-core --features op-contracts-strict`, exit 0; `/tmp/ck-valence-core-check-final.log`.
- Scoped runtime: `cargo test -p cosmolkit --no-default-features --release --features hydrogens,valence,rings,op-contracts-strict,cosmolkit-core/op-contracts-strict --test migration_h_remove_state_api --lib`, exit 0, 75 unit plus 7 public tests passed, 0 failed/ignored; `/tmp/ck-valence-focused.log`. This explicitly isolated diagnostic is NOT full runtime/workspace acceptance.
- Scoped runtime check: `cargo check -p cosmolkit --no-default-features --features hydrogens,valence,rings,op-contracts-strict`, exit 0; `/tmp/ck-valence-runtime-check-scoped.log`.
- Default-feature public/cache focused tests initially exited 101, blocked by concurrently unfinished DRAW registry references to `ops::depict` and `Coordinate2DParams`, not by a failing hydrogen assertion. Logs `/tmp/ck-valence-api.log`, `/tmp/ck-valence-cache.log`.
- Default-feature full runtime strict suite exited 101 with the same three missing DRAW symbol errors; `/tmp/ck-valence-runtime-full.log`.
- Required default-feature `cargo test -p cosmolkit --release --test migration_run_privacy` exited 101 at the same compilation boundary; `/tmp/ck-valence-privacy.log`. No privacy-suite pass is claimed.
- `cargo fmt --all` ran, followed by a successful `cargo fmt --all -- --check`; no unrelated semantics were intentionally edited.

The default runtime check also failed on those DRAW symbols. No placeholders,
feature changes, weakened privacy tests or repairs to the other agent's DRAW
work were made to hide these failures.

Full workspace command: `cargo test --workspace --release --features
cosmolkit/op-contracts-strict,cosmolkit-core/op-contracts-strict`. It emitted
the same three missing DRAW symbol errors for both runtime lib and lib-test,
then continued waiting for other build jobs. The supervisor interrupted this
own invocation with SIGINT; actual exit **130**, not normal completion or a
pass. Log `/tmp/ck-valence-workspace.log`. Full runtime/privacy/workspace gates
must be rerun once the concurrent DRAW integration compiles.
