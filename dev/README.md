# COSMolKit Development Operating Manual

This file is the entry point for development rules in `dev/`. Normative
standards and current design documents remain at this directory's top level.
Execution plans, reports, tools, and historical records live in named
subdirectories so that an old plan cannot be mistaken for current policy.

## Directory Map

| Path | Role |
|---|---|
| `dev/*.md` | Normative standards, current design decisions, and stable scope documents |
| [`plans/`](./plans/) | Plans that still contain executable, unchecked steps |
| [`gap_reports/`](./gap_reports/) | Audits, source inventories, validation reports, and blocker records |
| [`tools/`](./tools/) | Development-only audit, generation, and debugging tools |
| [`archive/`](./archive/) | Superseded roadmaps and inactive execution records; never normative |

## Normative Standards

- [`agent_plan_standard.md`](./agent_plan_standard.md): required structure for
  executable plans.
- [`policy_invariants.md`](./policy_invariants.md): public behavior and
  correctness invariants.
- [`source_reproduction_protocol.md`](./source_reproduction_protocol.md):
  source-port reproduction and marker rules.
- [`source_bisection_debugging_protocol.md`](./source_bisection_debugging_protocol.md):
  source-backed mismatch isolation protocol.
- [`operation_system_standard.md`](./operation_system_standard.md): operation
  registry, capability, COW, and state-migration standard.
- [`testing_contract.md`](./testing_contract.md): policy-to-test mapping.
- [`parity_testing_contract.md`](./parity_testing_contract.md): external
  reference parity-test contract.
- [`repository_organization_policy.md`](./repository_organization_policy.md):
  test and generated-data ownership rules.

## Current Design

- [`topology_operations.md`](./topology_operations.md)
- [`derived_effects_permission_model.md`](./derived_effects_permission_model.md)
- [`inplace_operation_api_design.md`](./inplace_operation_api_design.md)
- [`bio_structure_operation_contract_design.md`](./bio_structure_operation_contract_design.md)
- [`bio_structure_io_policy.md`](./bio_structure_io_policy.md)
- [`pdb_mmcif_gemmi_primary_plan.md`](./pdb_mmcif_gemmi_primary_plan.md)
- [`rdkit_pdb_molecule_conversion_plan.md`](./rdkit_pdb_molecule_conversion_plan.md)
- [`wwpdb_macromolecular_stress_experiment.md`](./wwpdb_macromolecular_stress_experiment.md)
- [`tetrahedral_stereo.md`](./tetrahedral_stereo.md)
- [`confseq_fast_geometry_design.md`](./confseq_fast_geometry_design.md)
- [`ai_native_features.md`](./ai_native_features.md)

Current parity boundaries are summarized in
[`VALIDATION.md`](../VALIDATION.md), and the current per-feature ledger is
[`porting_inventory.md`](./porting_inventory.md). Historical percentage
estimates and phase roadmaps are retained under
[`archive/roadmaps/`](./archive/roadmaps/) and are not current status claims.

## Required Development Mode

Development and CI must use strict operation checks:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --release --features op-contracts-strict
```

`op-contracts-strict` currently means:

```text
runtime-invariants + op-contracts
```

This enables:

```text
Molecule invariant checks
OpParts operation contract checks
operation body derived-effect and trace checks
```

Agents must not validate topology-operation changes with plain
`cargo check -p cosmolkit-core` alone. Plain checks are useful for a quick syntax
pass, but they are not sufficient for development signoff.

Use release test builds for routine full local validation. The strict feature
set still compiles runtime invariants and operation-contract checks; release
mode only avoids paying debug-build execution cost in parity tests where runtime
dominates compile time.

Small, focused tests may be run in the default debug profile while iterating on
one function or one module:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict <test-filter>
```

Large local test runs, parity suites, and CI test runs should use release mode
with the same strict feature set:

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict
```

Do not interpret release-mode testing as relaxed testing. Constraint coverage is
controlled by `op-contracts-strict`, not by the Cargo optimization profile.

Expensive exhaustive parity matrices are marked `#[ignore]` and run explicitly
in CI coverage jobs. Local default tests keep representative rows so everyday
validation remains bounded.

## Required Release Mode

Release builds should use the default feature set unless the human author
explicitly requests additional runtime checks:

```bash
cargo build --release
```

Release default intentionally omits internal contract checks:

```text
OpParts::finish does not run op-contract validation
Molecule::from_blocks does not run full molecule invariant scans
```

This is intentional. The invariant and contract system is a development and CI
guardrail, not a user-facing validation API.

## Feature Flag Policy

The current core checking features are:

| Feature | Meaning |
|---|---|
| `runtime-invariants` | Compile in molecule structural invariant checks |
| `op-contracts` | Compile in operation contract checks inside `OpParts` and `OpParts::finish` |
| `op-contracts-strict` | Development/CI mode; enables both checks |

Rules:

```text
Use op-contracts-strict for development and CI.
Use release default for published optimized builds.
Do not add runtime checking to release default without human approval.
Do not remove strict checks from CI to work around failures.
```

If strict mode fails, the correct response is to fix the operation contract,
registry metadata, or invariant logic. Agents are not allowed to weaken feature
flags, bypass `OpParts`, or suppress invariant failures without
confirming the design exception with the human author.

## Operation Registry Workflow

Topology-related public operations must be registered through `molecule_ops!`.

Every generated `MoleculeOpSpec` contains:

```text
method
impl_fn
domain
kind
topology_edit
access
may_mutate
auto_remap
derived_effects
semantic_preconditions
requires_mapping
support
parity
io_roundtrip
```

Registry entries explicitly provide `method`, `impl_fn`, `kind`, `access`,
`derived_effects`, `feature`, `parity`, and `invariant_profile`. Every parity
policy other than `not_applicable` also requires `parity_profile`. Other
inputs have defined macro defaults, although behaviorally meaningful
non-defaults should remain explicit for review.

For molecule operations, `derived_effects` has four independent,
pairwise-disjoint categories: `recompute`, `preserve`, `invalidate`, and the
narrowly allow-listed `operation_defined`. The last category delegates the
mechanism for a source-required transition, not its correctness, and is
currently restricted to hydrogen-removal valence state. Derived-cache read
authority still comes from block-level `access`; no effect category grants read
permission, and there is currently no per-cache-item read capability.
Unsupported source behavior is returned as a structured operation error rather
than encoded as a derived effect. The complete model is defined in
[`derived_effects_permission_model.md`](./derived_effects_permission_model.md).
Here, unsupported denotes an explicitly separate capability outside the
declared support boundary, never a mismatching case inside a parity-covered
surface.

The generated matrices are the source of truth:

```text
MOLECULE_OPS
SUPPORT_MATRIX
OPERATION_INVARIANT_MATRIX
PARITY_MATRIX
```

Do not maintain separate hand-written operation lists unless they are generated
from these matrices.

Registry coverage may grow incrementally while porting functionality. It does
not need to cover every future operation before porting starts. However, a
public operation must be registered before its real implementation is ported.
Do not add a new topology-related public method first and retrofit registry
metadata later.

## Operation Body Capability Object

Operation bodies receive a unified `OpParts` value through `#[mol_op_body(...)]`.
Do not add per-operation parts structs or manual capability combinations.

`OpParts` and its working molecule state must live in the private molecule-op
runtime module. Operation body implementation functions must live in sibling
modules, not in that runtime module. This is an architectural boundary, not a
layout preference: Rust privacy must make `parts.working` and every other
`OpParts` field unreachable from operation bodies.

`OpParts` is a fast operation capability object:

```text
OpParts::new clones Molecule cheaply by cloning Arc-backed blocks.
It must not eagerly clone topology, coordinates, properties, conformers, caches,
or other large molecule data.
Large block cloning happens only through registry-checked begin methods such as
`begin_topology_mut()`.
```

Operation bodies must access mutable molecule state only through scoped
`OpParts` methods:

```text
with_topology_mut(...)
with_coordinates_mut(...)
with_topology_and_properties_mut(...)
with_topology_coordinates_properties_mut(...)
record_topology_edit(...)
record_topology_mapping(...)
set_*_cache(...)
clear_cache(...)
prove_preserved(...)
```

Low-level `begin_*_mut()` / `commit_*()` methods implement and test the scoped
lifecycle; fallible operation bodies must not leave a block checked out across
an error. The macro-generated wrapper, not the operation body, owns `finish()`.

Operation bodies must not recover a raw `Molecule` or raw `&Molecule` through
helper APIs. Helpers called from operation bodies must accept narrowed inputs
such as `MoleculeReadParts`, `&[Atom]`, `&[Bond]`, `CoordinateBlock`, or an
explicit operation-local assignment/update plan. Do not add compatibility traits
or generic helper signatures that also accept `&Molecule`.

Whole-molecule algorithms used by operations must be split at the operation
boundary: compute from read-only narrowed inputs, then apply returned block
updates through `OpParts` begin/commit methods. Distance-geometry conformer
generation is the reference pattern: operation bodies call read-only coordinate
update helpers and write back only the returned `CoordinateBlock`.

`OpParts` is wrapper-owned state migration and contract-recording machinery. Do
not add chemistry, perception, sanitization, or operation-specific behavior to
`OpParts`. New operation behavior belongs in the operation `impl_fn` or in a
domain module called by that `impl_fn`; `OpParts` should only apply state
changes, perform COW access, remap registered blocks, record topology mapping
artifacts, record derived-state effects, and validate the operation contract.

`begin_topology_mut()` / `commit_topology(...)` are used for local
topology-state edits that do not change atom or bond table length, ordering, or
identity. Examples include changing bond order, aromatic flags, charges,
isotope state, query state, or other local atom/bond attributes.

The same boundary applies to `BioOpParts`: BioStructure operation bodies should
compute selection/transform intent and call explicit migration primitives such
as `remove_residues(...)`. Do not hand-write compacting hierarchy edits by
mutating atom, residue, chain, model, or coordinate rows directly from an
operation body.

## PDB/mmCIF And Protein IO Boundary

COSMolKit uses one public structural model: `BioStructure`.

Gemmi-derived code is the canonical source for PDB/mmCIF reading into
`BioStructure` and mmCIF writing from it. RDKit-derived PDB code is reserved for
`Molecule` compatibility behavior such as PDB block writing, `Molecule::from_pdb_block`,
`Molecule::from_mmcif_block`, and `AtomPDBResidueInfo` semantics. Molecule input
is layered over the Gemmi-primary `BioStructure` parse plus an explicit
structure-to-molecule conversion policy; it must not introduce a second public
parser.

Do not expose parallel public parser modules for the same user task. In
particular, do not publicize `pdb_parser.rs` or `mmcif_parser.rs` beside
`io::bio`; those files are quarantined until a source-backed Molecule-specific
need exists. The detailed rule is in
[`bio_structure_io_policy.md`](./bio_structure_io_policy.md) and the execution
plan is in
[`pdb_mmcif_gemmi_primary_plan.md`](./pdb_mmcif_gemmi_primary_plan.md).

Compacting, appending, or renumbering topology edits must use the operation
begin/commit discipline. Operation bodies begin each registry write-owned block
once, mutate the owned local blocks, commit them once, and record the declared
topology edit, topology mapping artifact, and derived-state effects
through `OpParts`. They must not mutate `Molecule` internals directly or mix
whole-molecule read views with independently writable blocks.

The registry contract for strong topology edits must declare the state migration
surface explicitly:

```text
topology_edit: compacting
auto_remap: [coordinates]
requires_mapping: required
```

In development and CI, `OpParts` checks `MoleculeOpSpec::may_mutate` with
assertions. Attempting to mutate a block not declared in the registry is a
developer/agent bug and should panic under `op-contracts-strict`.

In release default, permission assertions, trace recording, and contract
validation are not compiled. The release path still uses the same `OpParts`
accessors and COW storage model; it only omits development checks.
