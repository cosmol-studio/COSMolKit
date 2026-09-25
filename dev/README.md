# COSMolKit Development Manual

This is the entry point for development documentation. Follow
[AGENTS.md](../AGENTS.md) for agent authority and the rules linked below for
the work in scope. This index does not define a second architecture or queue.

## Authority and execution

- [Crate architecture](./crate_architecture.md): final ownership, dependency
  direction, and transaction/capability shape.
- [Public API design](./public_api_design.md): canonical names, receivers,
  value/in-place forms, and language projections.
- [Split-crate plan](./plans/crate_architecture_completion_plan.md): the sole
  execution queue and progress ledger. Reports and historical plans are
  evidence, not competing queues.
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

## Domain designs and protocols

- [BioStructure operation contracts](./bio_structure_operation_contract_design.md)
- [BioStructure IO policy](./bio_structure_io_policy.md)
- [wwPDB stress protocol](./wwpdb_macromolecular_stress_experiment.md)
- [Tetrahedral stereo](./tetrahedral_stereo.md)
- [MolAlign API design](./rdkit_molalign_api_design.md)
- [ConfSeq FastGeometry design](./confseq_fast_geometry_design.md)

Design documents describe their declared boundaries, not blanket implementation
or validation status. Source inventories and actual test results belong in unit
reports; acceptance belongs in the sole plan. Historical validation claims in
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
| [plans/](./plans/) | The sole split-crate queue and subordinate source-reference plans |
| [gap_reports/](./gap_reports/) | Audits, source inventories, validation and blocker evidence |
| [tools/](./tools/) | Development-only checking and preparation tools |
| [archive/](./archive/) | Historical snapshots and superseded records; never normative |

[VALIDATION.md](../VALIDATION.md) records the scope of reference evidence.
Historical inventories are not a second current-status ledger.
