# Crate Split Architecture Experiment

This is an isolated architecture experiment. It is not part of the root
workspace and does not depend on any published or local COSMolKit crate.

The experiment validates crate ownership and call direction only; the domain
functions intentionally contain minimal behavior and are not chemistry ports.

The runtime also reproduces the operation-contract shape used by the main
project: every exercised surface has a `MoleculeOpSpec`, support/invariant/
parity matrices, declared block access, mutation/remap/effect metadata,
semantic preconditions, mapping requirements, and a checked `OpParts` finish.
Multiple-output operations finalize each emitted branch independently before
returning the ordered `Vec<Molecule>`.

```text
runtime
  Molecule, private OpParts, macro-generated wrappers, composition
   |\
   | +--> batch ----> io / chemistry ----> model
   +----> io / chemistry -----------> model
```

The isolated `macros` crate supplies the same three compile-time surfaces used
by the production design: `molecule_ops!` generates operation specs, registries,
matrices, scalar wrappers, and in-place wrappers; `mol_op_body` injects the
single-output `OpParts` capability; and `mol_multi_op_body` injects the
multi-output capability.

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
receive `Molecule` or `OpParts`; the runtime extracts detached model values and
installs results only after the algorithm returns successfully.

The batch crate owns only ordered parallel scheduling over detached values. The
runtime `MoleculeBatch` owns the public execution policy (`n_jobs` and progress
display), and each parallel molecular transform still invokes the scalar
`Molecule` operation so every record passes through its own operation-contract
lifecycle.

Construction functions and ordinary block transformations return their
explicit model values directly as `(TopologyBlock, CoordinateBlock,
MoleculeProperties)`. Multiple-output algorithms return a vector of those
tuples; the runtime validates each branch and constructs the public `Molecule`
only at the trusted boundary. This experiment deliberately has no generic
draft container: such a type would imply a partially valid molecule
without adding a real contract or ownership guarantee.

When opening the repository root in Zed, the root `.zed/settings.json` links
this nested Cargo workspace to rust-analyzer. Restart the language server after
changing the linked-project configuration if the experiment types are not yet
indexed.
