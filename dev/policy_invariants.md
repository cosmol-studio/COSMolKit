# COSMolKit Development Principles

COSMolKit is a Rust-native cheminformatics toolkit built through source-backed
reproduction of established libraries such as RDKit and Gemmi. It aims to
combine faithful chemical behavior with explicit APIs, efficient execution,
and reliable Python, JavaScript/WASM and machine-learning workflows.

> Correctness first, explicit behavior second, performance third.

This is the project's order of priorities, not permission to ignore performance.
A result must be chemically correct, its meaning and limitations must be clear,
and its implementation must avoid unjustified costs. Neither a convenient API
nor a faster benchmark justifies silently changing behavior.

These principles apply to human contributors and agents alike. They describe
how to make development decisions, what shortcuts are unacceptable, and what
evidence is needed before calling work complete.

## Reproduce behavior; do not fit examples

A port starts from the pinned upstream implementation, not from an intuitive
description of the algorithm or a handful of matching outputs. Preserve its
control flow, state transitions, options, error handling and source-defined
fallbacks across the declared scope.

- Inspect the existing implementation before rewriting it. Reuse correct work;
  audit questionable work against the source rather than inheriting its claims.
- Do not replace a difficult branch with a heuristic, default value, special
  case for one molecule, or silent fallback.
- When results differ, locate the first semantic divergence. Do not patch the
  final output until the failing example happens to pass.
- Keep source evidence beside the implementing code. Review behavioral fidelity
  and performance/complexity separately; neither proves the other.
- An intentional difference requires explicit approval and a documented reason
  and scope. A surprising upstream behavior is not permission to redesign it.

The [source reproduction protocol](./source_reproduction_protocol.md) defines
source anchors and review markers. The
[source bisection protocol](./source_bisection_debugging_protocol.md) explains
how to locate a divergence without heuristic patches.

## Define the contract before implementing the feature

Decide what a capability means, which inputs and options it supports, what it
returns, and how it fails before implementation. Define planned interfaces in
the implementation plan; the public registry describes real APIs, not missing
functions. Deliver the implementation and its public registration together,
with any required operation contract declared before its body.

Public names must describe the user's task, not leak an internal helper's name
or the organization of an implementation crate. Language bindings project the
same behavior; they do not own alternative chemistry algorithms.

Upstream-derived chemistry is expected to match its pinned reference by
default; parity is not an optional premium level of support. Departures need
explicit approval and a precise explanation. Original functions follow their
own project-defined behavior. Experimental functions must identify their
limitations rather than imply a settled behavior contract.

A function's status describes its behavior commitment, not a test result or an
execution permission. Contributors update declarations explicitly after
reviewing evidence; tests do not automatically promote them. Registration and
compilation do not prove chemical correctness. Declaring Python/JS projections
does not establish that those bindings have been implemented.

The [public API design](./public_api_design.md) defines naming, receivers,
types and cross-language projections.

## Keep one owner for each responsibility

A local task must not make the project architecture less coherent. Algorithms
belong to their domain owners; the public runtime coordinates authorized work
and owns live molecule state. Bindings remain projections of the public API.

Reuse the shared mechanisms for lifecycle and state management. Do not copy
an algorithm into a wrapper, create a second model to bypass a boundary, or
grant an operation broader access simply because that makes implementation
easier. Shared guarantees belong in the framework, not in repeated handwritten
bookkeeping inside every algorithm.

An abstraction must serve a concrete responsibility. Do not add compatibility
layers, generic frameworks, duplicate registries or temporary adapters merely
to make a local step easier. Changes to ownership or authority require an
explicit design decision.

The [crate architecture](./crate_architecture.md) assigns ownership and
dependency direction; the [operation standard](./operation_system_standard.md)
defines operation contracts. The
[architecture rationale](./architecture_rationale.md) explains why these
boundaries matter, particularly for agent-assisted development.

## Make errors and information loss visible

An unavailable independent capability must fail explicitly. A failing input or
branch inside a supported capability is a defect, not a new unsupported
category. Never turn an error into an empty molecule, an empty collection or a
chemically plausible placeholder unless that result is genuinely the specified
behavior.

Errors retain their known category, underlying reason and available input
context. Batch failures remain traceable to their original input indices.
Failed items must not disappear unless the caller selects a documented skip
policy, and omissions must remain observable.

Batch execution preserves input order unless another order is explicitly
declared. Parallel scheduling must not change the specified results, ordering
or error positions for identical inputs and options, including the seed where
relevant.

Single-source value transformations preserve molecule-level metadata by
default. Splitting and multi-source operations must declare metadata ownership
and conflict handling. Each I/O path must state what it preserves and what the
format cannot represent; successful parsing or writing does not prove a
lossless roundtrip.

## Represent meaning explicitly

Chemical semantics belong in typed state, not magic string properties.
Supported enhanced stereo groups, SGroups, attachment points, brackets,
cstate vectors and V3000 collections must retain their modeled meaning.
Unknown non-semantic metadata may remain raw properties.

Persistent references use COSMolKit identifiers. Upstream row numbers,
bookmarks and sequence IDs are source metadata, not interchangeable local IDs.

SDF property-list expansion retains the raw field and represents interpreted
atom/bond values explicitly. Target counts must match after source-defined
empty-value handling. A mismatch produces a structured error in strict mode;
non-strict mode must preserve the raw field without expansion or explicitly
report an unsupported capability path. Silent row loss or assignment to the
wrong atom or bond is unacceptable.

Graph, tensor and other machine-consumed exports must document ordering,
dimensions, feature meanings, coordinate/chirality conventions, missing values
and schema version. Breaking schema changes require a new version.

## Treat performance as part of implementation quality

Source reproduction includes data structures, asymptotic complexity,
allocation behavior and hot-path costs. Matching values is insufficient when
the translation introduces avoidable repeated scans, whole-object clones or
extra buffering.

Do not optimize by weakening error handling, changing chemistry or making
correctness depend on development-only checks. Strict and ordinary builds
must execute the same chemistry. Performance claims need an explicit basis:
source-level cost analysis or measurements against a named baseline, not an
isolated elapsed time.

Unsafe code must have a necessary, localized purpose and documented soundness
assumptions. It is not an escape hatch for an inconvenient ownership boundary.

## Finish complete behaviors and report evidence honestly

A completed task delivers its entire specified behavior, not a dispatch stub,
a partial branch, copied source comments, or a TODO. If the task is too large,
split it into smaller complete behaviors with immediate validation; do not
silently reduce its scope.

Regression tests verify focused local behavior. Corpus parity tests compare
registered public operations with a pinned reference across selected inputs.
Both are evidence; neither replaces source review. Their placement and workflow
are defined in [test boundaries](./test_boundaries.md).

After changing tests, run the most specific relevant tests immediately. Final
integration validation does not replace that check. Preserve counterexamples,
and never obtain a pass by dropping cases, weakening comparisons or changing
expectations without source evidence.

Report actual commands, scope, counts, exit codes and remaining failures.
Compilation alone, zero matching tests, ignored suites and unrun checks are
not successful validation. A known failure remains a defect and blocks a
parity claim for the affected scope. Do not present a narrow gate as whole-project
acceptance.

## Work within the authorized scope

Read the applicable rules before acting. Review requests are read-only unless
changes are requested. Preserve other contributors' work; do not clean up
unrelated files or revise standards to accommodate an implementation.

For split-crate execution, the
[completion plan](./plans/crate_architecture_completion_plan.md) is the sole
queue and progress ledger. Follow its steps and required readings in order.
Record progress in the plan so it survives context compaction; do not rely on
conversation memory. Continue authorized work until completion, a genuine
blocker or interruption, rather than stopping after an arbitrary batch.

Do not ask for permission to use shortcuts already prohibited by project rules.
Resolve ordinary implementation and validation failures within scope. Stop
for a decision when proceeding would require new authority, an unresolved
design choice or a change to the agreed boundary.

Git operations, including read-only inspection, commits and pushes require
explicit authorization. Building does not authorize publishing or version
changes. Repository documentation is English; generated corpora, caches and
temporary diagnostics must not accidentally become committed source.

[AGENTS.md](../AGENTS.md) defines agent authority,
[the plan standard](./agent_plan_standard.md) defines auditable steps, and
[repository organization](./repository_organization_policy.md) defines artifact
placement. These rules keep progress reviewable without making local
implementation a license to redesign the project.
