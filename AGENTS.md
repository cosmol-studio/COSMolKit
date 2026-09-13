# COSMolKit Agent Guidelines

## Authoritative Documents

- Follow `dev/crate_architecture.md` for ownership, `dev/README.md` for development rules, and `dev/public_api_design.md` for public APIs.
- Read applicable rules before working. After every context compaction, reread `dev/source_reproduction_protocol.md`. Execute each plan task's preceding Read Step; prior reading is not a substitute.
- Use `dev/plans/crate_architecture_completion_plan.md` as the sole split-crate execution plan and ledger. Reports, experiments, and historical code are not alternative architectures or queues.
- Report conflicting rules with concrete evidence and ask for clarification; do not choose weaker rules. Code comments specifying operation constraints or human approval are binding.

## Agent Authority

- Reviews, diagnosis, and acceptance checks are read-only by default. Relevant non-destructive checks are allowed; do not change code, plans, or another agent's task without a request to do so.
- Implement only the authorized scope. Follow the plan sequentially until completion, a genuine blocker, or interruption; do not skip steps or stop merely after a time interval or batch.
- If a task is too large or lacks dependencies, split complete behaviors and immediate validation steps in the sole plan under `dev/agent_plan_standard.md`. Stubs, TODOs, and partial implementations are not completion.
- Preserve others' changes. Do not overwrite, revert, or clean unrelated files; coordinate overlapping edits.
- Without explicit authorization, do not perform Git operations (including read-only commands), commit/push, publish, bump versions, bulk-update dependencies, or perform destructive cleanup.
- Do not change AGENTS.md, architecture, execution standards, or weaken checkers to accommodate implementation. Rule changes require explicit authorization. Correct test expectations only with source evidence; retain coverage and prior failures.
- Stop and ask before changing ownership, adding crates, introducing compatibility layers, expanding scope, or violating approval comments. Fix ordinary build and test failures within scope rather than treating them as automatic stopping points.

## Crate Ownership

- `cosmolkit-types`: foundational vocabulary. `cosmolkit-model`: detached values and local structural validation.
- `cosmolkit-core`: source-backed foundational algorithms; no `Molecule` or runtime ownership.
- Domain crates own notation, IO, search, descriptors, fingerprints, bio, and other algorithms as assigned by the architecture. Production dependencies must be acyclic and must not depend back on `cosmolkit`.
- `cosmolkit`: the sole live `Molecule`/builder, private runtime, contracts/cache, authorized extraction and validated commit, and thin public APIs. No domain algorithm implementations.
- `cosmolkit-macros`: declaration parsing and capability, wrapper, registry, and matrix generation; no domain algorithms.
- `python/` and `wasm/`: language projections of the public API, with chemistry dependencies only on `cosmolkit`. Do not expose runtime internals or reimplement algorithms.
- Each behavior has one owner. Reuse correct implementations, implement missing behavior in its owner, and expose it through canonical APIs without duplication.

## Operations and Public APIs

- Before changing core, operation bodies, registries, or macro machinery, read `dev/README.md` and the applicable operation standards.
- Declare topology/coordinate operation contracts through `molecule_ops!` before implementing bodies. Generate `MOLECULE_OPS`, `SUPPORT_MATRIX`, `OPERATION_INVARIANT_MATRIX`, and `PARITY_MATRIX` from declarations only.
- The sole operation runtime boundary is `ops::runtime::{context,multiple,registry}`. Keep operation bodies outside this subtree, using only generated marker-specific capabilities. Rust visibility follows semantic modules, not physical file directories.
- Bodies must not access internal fields, unrestricted runtime read/write primitives, or constructor/finish/abort entrypoints. Do not add handwritten parallel capability lists or a second registry.
- Runtime owns block-level COW, mappings, effects, invariants, and atomic commit. Domain crates exchange explicit detached values, never `Molecule`, `OpParts`, or commit authority.
- Preserve sharing of unchanged blocks and failure atomicity. Do not pass tests by cloning all state, swallowing errors, or bypassing validation.
- Register public capabilities in `crates/cosmolkit/src/binding_contract/` before exposure. Keep signatures, defaults, features, errors, and language projections consistent. Registration does not imply implementation or supported status.
- Non-mutating transforms return new values. In-place `Molecule` operations must end in `_`; that suffix has no other meaning. Queries return read results.
- Keep helpers private by default and cross-crate interfaces narrow. Do not expose mutable molecule storage or compatibility aliases.

## Source Reproduction

- Follow `dev/source_reproduction_protocol.md` against pinned upstream sources. Place verbatim source anchors inside the implementing function, with separate behavior and complexity markers.
- Do not hide missing behavior behind heuristics, defaults, or silent fallbacks. Reproduce only source-defined fallbacks and propagate other errors structurally.
- Mark unmodeled independent capabilities unsupported; never relabel failures on supported inputs as unsupported. Source comments or missing oracles are not completion evidence.

## Validation and Handoff

- After changing tests, immediately run the most specific relevant tests. Retain commands, exit codes, counts, and failures. When the plan specifies `run_unit_validation.py`, freeze targets/features first and do not bypass the runner.
- For core changes, run `cargo check -p cosmolkit-core --features op-contracts-strict` and `cargo test -p cosmolkit-core --release --features op-contracts-strict`.
- For runtime, operation, or macro changes, also run affected-crate tests, `cargo check -p cosmolkit --features op-contracts-strict`, and `cargo test -p cosmolkit --release --features op-contracts-strict`. Explicitly enable affected capability features so cfg gates do not hide coverage.
- For operation integration or visibility changes, run `cargo test -p cosmolkit --release --test migration_run_privacy`. Retain real-module-layout default/strict compile-pass and compile-fail cases; text checks or runtime rejection cannot replace compile-time isolation.
- Final cross-crate validation: `cargo test --workspace --release --features cosmolkit/op-contracts-strict,cosmolkit-core/op-contracts-strict`. Enable runtime strict on `cosmolkit`; the core feature is not a substitute.
- Use stage-specific exclusions only as prescribed by the plan. Do not report a stage pass as a full workspace pass or silently exclude members from final validation.
- Small focused debugging may use debug builds. Full, large, and parity suites use release with strict checks. Run `cargo fmt --all` after Rust edits.
- Release builds use default features unless additional checks are explicitly requested. Building does not authorize publishing.
- Report actual changes, validation, and remaining blockers. Unrun tests, zero matches, or entirely ignored suites are not passes. When authorized to commit, use Conventional Commits and honor file exclusions.

## Files and Tooling

- Use Rust 2024, four-space indentation, and standard naming. Keep public APIs in `lib.rs` or narrow modules.
- Put tests in the owning crate's unit modules or `tests/`, shared fixtures in `testdata/`, and preparation entrypoints in `tools/testdata/`. Follow `dev/repository_organization_policy.md`; do not casually commit generated data or bulk corpora.
- Manage Python from the repository root with `uv sync --group dev`. Build and test with `.venv/bin/maturin develop --manifest-path python/Cargo.toml` and `.venv/bin/pytest`.
- Generate stubs with `cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen`; do not hand-edit `python/cosmolkit.pyi`.
- Build docs with `.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html`; type-check with `.venv/bin/basedpyright python/tests python/examples`.
