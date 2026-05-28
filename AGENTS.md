# Repository Guidelines

Read `dev/source_reproduction_protocol.md` after every context compaction.
Before editing `crates/cosmolkit-core/`, the operation registry, topology
operations, or macro-controlled operation machinery, read `dev/README.md`.

## Project Structure

- `crates/cosmolkit-core/`: current molecular graph, state, operation, IO, and
  chemistry core.
- `crates/cosmolkit-core-old/`: previous implementation; reference material
  only. Do not depend on, re-export from, or silently copy its behavior into the
  current core.
- `crates/cosmolkit-macros/`: proc macros for operation bodies, registries, and
  generated tables.
- `crates/cosmolkit/`: Rust facade and top-level public API.
- `python/`: PyO3 bindings, Python packaging, examples, and docs.
- `tests/`: integration, regression, parity, and golden-test data.
- `dev/`: normative design documents and development operating rules.

Keep public APIs in `lib.rs` or narrow public modules. Keep implementation
helpers private unless there is a deliberate public API reason.

## Source Reproduction

Ports from RDKit, Gemmi, OpenBabel, or other C/C++ libraries must follow
`dev/source_reproduction_protocol.md`.

- Do not implement ported chemistry or structural-biology behavior by
  heuristic approximation.
- Copy the relevant original source lines into the corresponding Rust function
  as review anchors.
- Use the required two-axis source markers for behavior and complexity status.
- Unsupported source behavior must return structured unsupported errors or stay
  marked unsupported; do not hide gaps behind silent fallbacks.
- `crates/cosmolkit-core-old/` may explain prior intent, but it is not an
  authority for current behavior.

## Operation System

Topology-related public molecule operations must go through `molecule_ops!` and
the operation-contract machinery. Do not bypass `TopologyOnlyParts`, direct
mutable storage, generated registries, or invariant checks.

The generated registry and matrices are the source of truth:

```text
MOLECULE_OPS
SUPPORT_MATRIX
OPERATION_INVARIANT_MATRIX
PARITY_MATRIX
```

Mutation of public `Molecule` values must be explicit:

- Non-mutating methods return new molecule values.
- Public in-place `Molecule` operations must end with trailing `_`.
- No public `Molecule` method may use trailing `_` for any other meaning.

Code comments that define operation requirements, agent guardrails, or
human-author approval requirements are binding project rules. If a requested
change appears to require violating such a comment, stop and ask for explicit
human-author approval.

## Required Validation

For core, topology-operation, or operation-registry changes, use strict checks:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict
```

For cross-crate changes, also run:

```bash
cargo test --workspace --features cosmolkit-core/op-contracts-strict
```

Plain `cargo check -p cosmolkit-core` is acceptable only as a quick syntax pass.
Run formatting after Rust edits:

```bash
cargo fmt --all
```

Release builds should use default features unless runtime checks are explicitly
requested:

```bash
cargo build --release
```

## Python Tooling

Use project-level `uv` environment management from the repository root:

```bash
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
.venv/bin/pytest
```

Generate stubs with:

```bash
cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen
```

The generated file is `python/cosmolkit.pyi`; do not edit it by hand.

Build docs and run type checks with:

```bash
.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html
.venv/bin/basedpyright python/tests python/examples
```

Python APIs should follow the current Rust core API, should not expose internal
helpers, and should not preserve old-core compatibility behavior unless there is
an explicit public API decision.

## Coding And Testing

Use Rust 2024 defaults: 4-space indentation, `snake_case` functions/modules,
`CamelCase` types/traits, and `SCREAMING_SNAKE_CASE` constants.

Prefer narrow modules and explicit types. Do not expose mutable molecule storage
through public APIs.

Place unit tests near code with `mod tests`; use integration tests under
`tests/` when a cross-module fixture is needed. Name tests by behavior, for
example `kekulize_handles_fused_aromatics`.

For operation changes, tests should run under `op-contracts-strict`. Unsupported
behavior should return structured unsupported errors, not panic, silently no-op,
or emit placeholder chemistry.

## Commits

Use Conventional Commits, for example:

```text
feat: add sdf streaming parser
fix: correct aromatic valence handling
```

Keep commits focused and buildable where practical. PR descriptions should
include scope, validation steps, fixture notes, linked issues when applicable,
and support-status or parity-policy changes.
