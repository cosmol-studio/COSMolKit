# COSMolKit Development Operating Manual

This file is the entry point for development rules in `dev/`.

It records operational details that agents and human contributors must follow
while editing the redesigned core. The design documents define principles; this
file defines the default commands and feature flags that should be used during
daily development.

## Related Documents

- [porting_plan.md](porting_plan.md) — source porting roadmap and current phase
  (entry point for porting work)
- [porting_inventory.md](porting_inventory.md) — per-feature porting status
- [work.md](work.md) — active phase implementation detail
- [bio_structure_io_policy.md](bio_structure_io_policy.md) — PDB/mmCIF
  structural IO boundary between Gemmi-derived BioStructure readers and
  RDKit-derived Molecule compatibility behavior
- [pdb_mmcif_gemmi_primary_plan.md](pdb_mmcif_gemmi_primary_plan.md) —
  execution plan for Gemmi-primary PDB/mmCIF reading and RDKit-compatible
  molecule construction
- [confseq_fast_geometry_design.md](confseq_fast_geometry_design.md) —
  ConfSeq FastGeometry backend boundary, target architecture, metrics, and
  validation rules

## Required Development Mode

Development and CI must use strict operation checks:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict
```

`op-contracts-strict` currently means:

```text
runtime-invariants + op-contracts
```

This enables:

```text
Molecule invariant checks
OpParts operation contract checks
operation body invalidation/handled-state checks
```

Agents must not validate topology-operation changes with plain
`cargo check -p cosmolkit-core` alone. Plain checks are useful for a quick syntax
pass, but they are not sufficient for development signoff.

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

Each registered operation must declare:

```text
method
impl_fn
domain
kind
may_mutate
derived_effects
requires_mapping
allows_noop
feature
parity
io_roundtrip
invariant_profile
parity_profile when parity is not not_applicable
```

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

`OpParts` is a fast operation capability object:

```text
OpParts::new clones Molecule cheaply by cloning Arc-backed blocks.
It must not eagerly clone topology, coordinates, properties, conformers, caches,
or other large molecule data.
Large block cloning happens only through registry-checked begin methods such as
`begin_topology_mut()`.
```

Operation bodies must access mutable molecule state only through `OpParts`
methods:

```text
begin_topology_mut()
commit_topology(...)
begin_coordinates_mut()
commit_coordinates(...)
begin_properties_mut()
commit_properties(...)
record_topology_edit(...)
record_topology_mapping(...)
recompute_aromaticity()
recompute_valence()
recompute_stereo()
clear_cache(...)
finish(...)
```

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
`BioStructure`. RDKit-derived PDB code is reserved for `Molecule` compatibility
behavior such as PDB block writing, future `Molecule::from_pdb_block` support,
and `AtomPDBResidueInfo` semantics. Future molecule input must be layered over
the Gemmi-primary `BioStructure` parse plus an explicit conversion profile; it
must not introduce a second public parser.

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
topology edit, topology mapping artifact, and cache/handled-state effects
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

## Related Documents

- [`policy_invariants.md`](./policy_invariants.md): project policy invariants
- [`testing_contract.md`](./testing_contract.md): policy-to-test mapping
- [`topology_operations.md`](./topology_operations.md): topology operation rules
- [`operation_system_standard.md`](./operation_system_standard.md): normative
  operation registry, `OpParts`, COW, and state-migration standard
- [`PARITY_TESTING_CONTRACT.md`](./PARITY_TESTING_CONTRACT.md): RDKit parity testing contract
