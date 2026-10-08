# COSMolKit Development Manual

This is the entry point for development documentation. Follow
[AGENTS.md](../AGENTS.md) for agent authority and the rules linked below for
the work in scope. This index does not define a second architecture or queue.

## Authority and execution

- [Crate architecture](./crate_architecture.md): final ownership, dependency
  direction, and transaction/capability shape.
- [Public API design](./public_api_design.md): canonical names, receivers,
  value/in-place forms, and language projections.
- [Task plans](./plans/README.md): execute only the current user-authorized
  task and its explicitly assigned plan. Dated receipts are evidence, not
  blanket validation of the current checkout or an automatic work queue.
- [Agent plan standard](./agent_plan_standard.md): numbered Read/action steps,
  immediate validation, and evidence requirements.
- [Architecture rationale](./architecture_rationale.md): non-normative design
  explanation, not an additional per-step reading requirement.

An old path, code example, or historical completion claim cannot override the
architecture. Preserve source evidence and approved behavior when reconciling
documentation; do not weaken checks to accommodate an implementation.

## Normative standards

| Document | Responsibility |
|---|---|
| [Policy invariants](./policy_invariants.md) | Observable behavior and correctness promises |
| [Operation standard](./operation_system_standard.md) | Registry, generated capabilities, COW, mappings, commit and failure semantics |
| [Derived effects](./derived_effects_permission_model.md) | Cache effects and their separation from read authority |
| [Source reproduction](./source_reproduction_protocol.md) | Pinned-source anchors and independent behavior/complexity review |
| [Source bisection](./source_bisection_debugging_protocol.md) | First-divergence debugging without heuristic patches |
| [Test boundaries](./test_boundaries.md) | Fixed regressions, corpus parity, and concrete test methods |
| [Repository organization](./repository_organization_policy.md) | Tests, fixtures, generated data, and preparation tooling |

Operation work follows the [strict/release validation requirements](./operation_system_standard.md#19-strict-and-release-builds).
Enable runtime strict checks on `cosmolkit`; core strict alone is not a runtime
gate. Release optimization and strict checking are independent. Published
builds use default features unless extra checks are explicitly requested;
building does not authorize publishing.

## Pre-commit checks

Run from the repository root against the final changes. The command below
excludes both parity packages, including `reference_parity`; run corpus and
special-regression suites separately through their prepare/test entrypoints.

```bash
uv sync --locked --group dev
cargo fmt --all --check
cargo test --workspace --locked --profile dev-test --no-fail-fast \
    --exclude cosmolkit-parity-tests --exclude cosmolkit-parity-tests-fixed \
    --features cosmolkit/op-contracts-strict
cargo run -p cosmolkit-py --no-default-features --features dev-stub --bin stub_gen
.venv/bin/maturin develop --profile dev-test --manifest-path python/Cargo.toml
python3 wasm/tools/wasm_binding/run.py
.venv/bin/pytest python/tests
```

Also run every `run` step in the `feature-matrix` job of
[features.yml](../.github/workflows/features.yml); do not maintain a separate
feature list or weaken CI.

Daily `--profile dev-test` builds use optimization level 3 without LTO and with
16 codegen units. CI distribution builds use `--release` for fat LTO and one
codegen unit; strict checks are independent of either profile.

All checks must pass before committing. Record commands, exit codes, test
counts and logs; do not add exclusions beyond the parity-package separation
above, filters or new skips to hide failures.
Test the freshly built Python extension, not an older installation. Rerun
affected checks after further changes. This checklist does not authorize Git
operations.

## Domain designs and protocols

- [Double formatting](./double_formatting_contract.md): approved Boost-compatible pure-Rust binary64 string conversion and its validation boundary.

- [Coordinate storage and selection](./coordinate_selection_contract.md): approved separate 2D/3D storage and unique-or-explicit coordinate selection.
- [BIO architecture and lightweight operations](./bio_architecture.md): public BioStructure/Protein, detached data, IO ownership and generated COW operations.
- [wwPDB stress protocol](./wwpdb_macromolecular_stress_experiment.md)
- [Tetrahedral stereo](./tetrahedral_stereo.md)
- [MolAlign API design](./rdkit_molalign_api_design.md)
- [ConfSeq FastGeometry design](./confseq_fast_geometry_design.md)

Design documents describe their declared boundaries, not blanket implementation
or validation status. Source inventories and actual test results belong in unit
reports; acceptance belongs in the authorized task's plan. Historical validation claims in
domain documents must be read in their original scope, not transferred to a
new implementation.

## Deferred proposals

- [Coordinate model improvements](./coordinate_model_improvement_draft.md):
  explicitly deferred pending final review and approval.
- [AI-native features](./ai_native_features.md): future capability sketch.

Neither proposal authorizes implementation or reorders the active plan.

## Directory map

| Path | Role |
|---|---|
| `dev/*.md` | Canonical standards, domain designs, and clearly labeled rationale/proposals |
| [plans/](./plans/) | Task-specific plans and dated source-port records; not an automatic execution queue |
| [gap_reports/](./gap_reports/) | Audits, source inventories, validation and blocker evidence |
| [tools/](./tools/) | Development-only checking and preparation tools |
| [archive/](./archive/) | Historical snapshots and superseded records; never normative |

[VALIDATION.md](../VALIDATION.md) records the scope of reference evidence.
Historical inventories are not a second current-status ledger.
