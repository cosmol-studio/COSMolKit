# Crate Split Architecture Experiment

This is an isolated architecture experiment. It is not part of the root
workspace and does not depend on any published or local COSMolKit crate.

The experiment validates crate ownership and call direction only; the domain
functions intentionally contain minimal behavior and are not chemistry ports.

Its reduced model is not a reduced production contract. The model currently
checks only local bond-endpoint ranges because the experiment omits the full
atom, coordinate, property, adjacency, and cache data structures. Any real
crate migration must preserve the existing invariant checks at their current
strength: atom/bond identity and index validity, topology references and
adjacency consistency, conformer coordinate-row alignment, per-atom and
per-bond property alignment, cache validity, and all operation-specific
preconditions. Detached blocks and new crate boundaries are not reasons to
drop or weaken those checks.

The runtime also reproduces the operation-contract shape used by the main
project: every exercised surface has a `MoleculeOpSpec`, support/invariant/
parity matrices, declared block access, mutation/remap/effect metadata,
semantic preconditions, topology-transition metadata, and a checked `OpParts`
finish.
Multiple-output operations finalize each emitted branch independently before
returning the ordered `Vec<Molecule>`.

```text
runtime (cosmolkit)
  Molecule, private OpParts, registry, public wrappers
       |       |       |       |
       v       v       v       v
      core    io    tautomer stereoisomers conformer
       |              |       |      |
       |              \       |      /
       |               \      |     /
                 model

batch (scheduler only) ----------------> model
```

The graph is intentionally a set of independent domain edges. The runtime
uses optional dependencies whose `dep:` activation is owned by the matching
fine-grained feature; enabling `tautomer` does not enable hydrogen, conformer,
or batch APIs. The batch crate has no parser or chemistry dependency: runtime
injects the parser callback, so `batch` remains an execution policy rather than
an alias for every other domain. `core` is the shared chemistry implementation
layer: it owns tightly coupled foundational behavior such as valence,
sanitization, and hydrogen transforms, while still accepting only detached
model values.

Fine-grained features and user-facing bundles are separate:

| Fine-grained API feature | Public surface | User bundle |
| --- | --- | --- |
| `smiles`, `smarts`, `sdf`, `pdb`, `xyz` | Individual construction/parser APIs | `common_api` |
| `hydrogens`, `tautomer`, `stereoisomers` | Individual chemistry operations | `chemistry_api` |
| `conformer` | Conformer operation | `3d_api` |
| `batch` plus `smiles` | Batch scheduling and SMILES construction | `batch_api` |

`full` is only the explicit composition of these bundles. A bundle is a user
convenience and must not be used as an implementation dependency that exposes
unrelated methods. The registry remains complete metadata for contract review;
feature gates control generated public wrappers and operation bodies, while
the corresponding algorithm dependency is compiled only when its fine feature
is selected. Hydrogen operations are exposed by runtime feature gates, but
their implementation is supplied by `core`; this is an implementation
dependency and does not make the core layer a public runtime API. `2d_api` is
currently an explicit empty reservation because this experiment does not
include the depiction domain.

The isolated `macros` crate supplies the same three compile-time surfaces used
by the production design: `molecule_ops!` generates operation specs, registries,
matrices, scalar wrappers, in-place wrappers, and operation-specific capability
types; `mol_op_body` and `mol_multi_op_body` inject a direct
`OpParts<'_, OperationAccess>` or `MultiOutputOpParts<'_, OperationAccess>`
parameter. Each capability is a zero-sized marker used as the `Access`
parameter of the runtime transaction object, so the generated impl exposes only
the block access declared by that operation. There is no wrapper context or
second transaction container.

The molecule-operation registry contains the five stateful or derived
`Molecule` operations exercised here: hydrogen addition/removal, tautomer and
stereoisomer enumeration, and conformer generation. SMILES/SMARTS and molecular
file readers remain construction APIs; they are intentionally not presented as
topology mutations merely to place them in the operation registry.

The exercised surface includes hydrogen transforms, SMILES/SMARTS parsing,
SDF/PDB/XYZ readers, tautomer and stereoisomer enumeration, conformer
generation, and a batch API. Run it from this directory:

```bash
cargo test --workspace
cargo tree --workspace
```

The important compile-time property is that the algorithm crates do not
receive `Molecule` or an unrestricted `OpParts`; operation bodies receive the
same transaction object parameterized by a generated access marker, whose
methods match the registry access declaration. The runtime still performs
lifecycle, final block-consistency, contract, and invariant checks; topology
mapping remains an algorithm-internal aid rather than an operation-body proof
registration. Compile-time visibility is an additional boundary, not a
replacement for those checks.

The batch crate owns only ordered parallel scheduling over detached values. The
runtime `MoleculeBatch` owns the public execution policy (`n_jobs` and progress
display), injects the selected parser, and each parallel molecular transform
still invokes the scalar `Molecule` operation so every record passes through
its own operation-contract lifecycle.

Construction functions and ordinary block transformations return their
explicit model values directly as `(TopologyBlock, CoordinateBlock,
MoleculeProperties)`. Multiple-output algorithms return a vector of those
tuples and call `MultiOutputOpParts::emit_all`; the runtime validates every
candidate and constructs the public `Molecule` only at the trusted boundary.
The multiple-output runtime deliberately has no branch-id, source-derivation,
or draft container because those would duplicate an algorithm’s own candidate
enumeration without adding an authority boundary.

When opening the repository root in Zed, the root `.zed/settings.json` links
this nested Cargo workspace to rust-analyzer. Restart the language server after
changing the linked-project configuration if the experiment types are not yet
indexed.
