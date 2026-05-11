# Repository Guidelines

## Current Project State
COSMolKit is in an intentional core redesign.

The redesigned implementation lives in `crates/cosmolkit-core/`. The previous
implementation has been moved to `crates/cosmolkit-core-old/` and is reference
material only. Do not depend on, re-export from, or silently copy old-core
behavior into the redesigned core.

Before editing the redesigned core, operation registry, topology operations, or
macro-controlled operation machinery, read `dev/README.md`.

## Project Structure
- `crates/cosmolkit-core/`: redesigned molecular graph, state, operation, IO,
  and chemistry core
- `crates/cosmolkit-core-old/`: old implementation, reference only
- `crates/cosmolkit-macros/`: proc macros for operation bodies, registries, and
  generated tables
- `crates/cosmolkit/`: facade re-exports and public top-level Rust API
- `python/`: PyO3 bindings and packaging metadata
- `tests/`: integration data, regression data, parity data, and test scripts
- `dev/`: normative design documents and development operating rules

Keep public APIs in `lib.rs` or narrow public modules. Keep implementation
helpers private unless there is a deliberate public API reason.

## Required Commands
For redesigned core work, development and CI validation must use strict
operation checks:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict
```

Plain `cargo check -p cosmolkit-core` is allowed as a quick syntax pass, but it
is not sufficient as final validation for core or topology-operation changes.

Release builds should use default features unless the human author explicitly
requests runtime checks:

```bash
cargo build --release
```

Run formatting after Rust edits:

```bash
cargo fmt --all
```

Workspace-wide `cargo check` may currently fail because the Python binding still
contains old-core API assumptions. Do not treat those failures as permission to
change core design back toward old parity behavior.

## Python Tooling
Use project-level `uv` environment management from repo root:

```bash
uv sync --group dev
```

Install the local extension into the shared environment:

```bash
.venv/bin/maturin develop --manifest-path python/Cargo.toml
```

Run Python-facing tests:

```bash
.venv/bin/pytest
```

Generate `.pyi` stubs with the dev stub feature:

```bash
cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen
```

The generated file is `python/cosmolkit.pyi`; do not edit it by hand.

Build a release wheel with the abi3-py39 feature:

```bash
.venv/bin/maturin build --release --manifest-path python/Cargo.toml --features release-abi3-py39 --out python/dist
```

Build Python documentation:

```bash
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
.venv/bin/python -m sphinx -b html python/docs/source python/docs/build/html
```

The generated documentation entry point is
`python/docs/build/html/index.html`.

Run Python type checking:

```bash
.venv/bin/basedpyright python/tests python/examples
```

The Python layer may lag behind the redesigned Rust core during bootstrap.
Fixing Python should follow the new core API, not restore old mutable internals.

## Coding Style
Use Rust 2024 defaults:

- 4-space indentation
- `snake_case` for functions and modules
- `CamelCase` for types and traits
- `SCREAMING_SNAKE_CASE` for constants

Prefer narrow modules and explicit types. Do not expose mutable molecule storage
through public APIs. Mutation of existing molecules must go through registered
operations and their operation-contract machinery.

## Code Guardrails
Code comments that define operation requirements, agent guardrails, or
human-author approval requirements are binding project rules. Do not bypass,
delete, weaken, or work around these comments while implementing a change.

If a change appears to require violating a comment such as "do not replace this
with raw indices in public APIs without human-author approval", stop and ask the
human author to confirm the design exception before editing that code.

## Operation Registry Rules
Topology-related public operations must be registered through `molecule_ops!`.
The generated registry and matrices are the source of truth:

```text
MOLECULE_OPS
SUPPORT_MATRIX
OPERATION_INVARIANT_MATRIX
PARITY_MATRIX
```

Do not maintain a separate hand-written operation list unless it is generated
from those matrices. Do not bypass `TopologyOnlyParts` for operation bodies.

## Testing Guidelines
Place unit tests near code with `mod tests` and integration tests under
`tests/` when a cross-module fixture is needed. Name tests by behavior, for
example `kekulize_handles_fused_aromatics`.

For core operation changes, tests should run under `op-contracts-strict`.
Unsupported behavior should return structured unsupported errors, not panic,
silently no-op, or emit placeholder chemistry.

RDKit parity policy is being redefined in `dev/`; do not add new hard-coded
parity rules to this file.

## Commit & Pull Request Guidelines
Use Conventional Commits when committing, for example:

```text
feat: add sdf streaming parser
fix: correct aromatic valence handling
```

Keep commits focused and buildable where practical.

PRs should include scope summary, validation steps, fixture notes, linked issues
when available, and notes on support-status or parity-policy changes.
