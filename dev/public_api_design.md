# Public API Design

This document is normative for the public COSMolKit API. It complements
[`crate_architecture.md`](./crate_architecture.md) and defines how the single
runtime-owned `Molecule` surface is named, classified, and exposed to Rust,
Python, and JavaScript/WASM users.

The goal is one logical chemistry API with idiomatic spelling in each target
language:

```text
Rust:       mol.molecular_weight()?
Python:     mol.molecular_weight()
JavaScript: mol.molecularWeight()
```

The spelling may follow language convention, but the operation, inputs,
outputs, defaults, behavior commitment, and error category must remain the same.

## 1. Public Boundary

`cosmolkit` is the only supported user-facing Rust crate and the only crate
that owns or accepts a live `Molecule`.

Language adapters must depend on `cosmolkit`, not `cosmolkit-core` or another
implementation crate. Domain crates operate on detached
`cosmolkit-model` values and are never part of the user-facing chemistry API.

The following are internal or implementation boundaries and must not appear in
the public language APIs:

```text
OpParts
operation capability markers
operation registry internals
derived-cache authority
runtime working-state containers
algorithm-crate Molecule parameters
```

There is exactly one authoritative `Molecule` type. A binding wrapper may own
an internal Rust `Molecule`, but it must not define a second chemistry model or
duplicate operation semantics.

## 2. API Categories

Every public item must belong to one of these categories.

| Category | Canonical location | Examples |
|---|---|---|
| Molecule construction | `Molecule` associated function | `from_smiles`, `from_sdf`, `from_inchi` |
| Molecule query | `&self` method | `num_atoms`, `molecular_weight`, `tpsa` |
| Value transformation | `&self` method returning a new value | `with_hydrogens`, `sanitize` |
| Explicit in-place transformation | Rust `&mut self` method with trailing `_` | `add_hydrogens_` |
| Serialization | `&self` method | `to_smiles`, `to_inchi`, `to_sdf` |
| Multiple-output operation | `&self` method returning a collection/result | `enumerate_tautomers`, `generate_conformers` |
| Query construction | domain module function | `search::parse_smarts` |
| Cross-molecule operation | domain function or named operator | `maximum_common_substructure`, `align_to` |
| Global vocabulary/metadata | module function or associated value | `version`, `element_info` |
| Parameters/options | dedicated public value type | `SmilesWriteParams`, `EmbedParams` |
| Results/reports | dedicated public value type | `Fingerprint`, `TopologyMapping`, `ValidationReport` |
| Errors | dedicated public error type | `SmilesError`, `DescriptorError`, `OperationError` |

An item must not be placed at module scope merely because its implementation
currently lives in a domain crate. Ownership of the public entry point follows
the semantic receiver, not the implementation file.

## 3. Receiver Rule

If an operation has one primary molecule input, it is a `Molecule` method.

```rust
let mol = Molecule::from_smiles("CCO")?;
let weight = mol.molecular_weight()?;
let smiles = mol.to_smiles()?;
let expanded = mol.with_hydrogens()?;
let coverage = mol.mmff_has_all_molecule_params()?;
```

The following module-level shape is not allowed for new public APIs:

```rust
calc_mol_wt(&mol, false)
mmff_has_all_molecule_params(&mol)
mol_to_smiles(&mol, &params)
```

The implementation may remain a free function inside an algorithm crate, but
the public `cosmolkit` method is the only supported entry point.

An operation may remain a module-level function when it has no natural
receiver, for example:

```rust
search::parse_smarts(text, &params) -> QueryGraph
maximum_common_substructure(inputs, &params) -> McsResult
version() -> &str
```

Functions that accept one molecule plus an independent reference should use a
named method when the receiver is clear:

```rust
mol.align_to(&reference, &params)
```

## 4. Naming Rules

Rust is the canonical naming source. Python preserves `snake_case`; JavaScript
and TypeScript adapters convert the same identifier to `camelCase` without
changing its semantic name.

### 4.1 Constructors

Use `from_*` for constructors that create a `Molecule` or another owned value:

```text
Molecule::from_smiles
Molecule::from_sdf
Molecule::from_inchi
```

Do not add `mol_from_*` names to the public `Molecule` API.

### 4.2 Serialization

Use `to_*` for conversion from a molecule to a serialized representation:

```text
mol.to_smiles()
mol.to_inchi()
mol.to_sdf(&params)
```

`mol_to_*` is an internal source-port name only. It must not be
the name of a new public method.

### 4.3 Queries and descriptors

Use the domain term directly, without `calc_`, `mol_`, or redundant `get_`:

```text
molecular_weight
exact_molecular_weight
molecular_formula
tpsa
num_rotatable_bonds
mmff_has_all_molecule_params
```

Boolean methods use `is_`, `has_`, or a domain-specific predicate. Methods may
return `Result<bool, Error>` when source behavior can fail; they remain methods
and must not be silently converted into properties in another language.

Coordinate-reading methods on `Molecule` are dimension-specific:

```rust
coordinates_2d(&self) -> Option<&[[f64; 2]]>
conformers_3d(&self) -> &[Conformer3D]
```

`coordinates_2d` borrows the first stored 2D conformer's atom-ordered rows;
absence is `None`, while a stored empty conformer is `Some(&[])`. It never
generates coordinates or falls back to 3D. `conformers_3d` borrows all stored
3D conformers in order, preserving their IDs and metadata. Neither accessor
permits mutation. Do not expose a public `Molecule::coordinates` block getter
or a `Molecule::conformers` mixed-dimension tuple in place of these methods.
The runtime's complete coordinate-block access remains private.

### 4.4 Transformations

Value-style transformations use `with_*` and return a new molecule or result:

```text
with_hydrogens
without_hydrogens
with_kekulized_bonds
with_3d_conformer
```

Rust in-place operations use the existing mandatory trailing underscore:

```text
add_hydrogens_
remove_hydrogens_
```

The trailing underscore must never be used for a non-mutating or unrelated
meaning. New APIs should prefer value-style methods unless in-place mutation is
required for a documented performance or compatibility reason.

Eligible value and in-place entry points share the same registered operation
implementation; they are not separate algorithm families. The generated
in-place name defaults to `{method}_`, with explicit semantic overrides such
as `with_hydrogens` -> `add_hydrogens_` and `without_hydrogens` ->
`remove_hydrogens_`. These names grant no mutable-storage access.
In-place COW and failure guarantees are defined in
[the operation standard](./operation_system_standard.md#in-place-execution-and-failure-semantics).

### 4.5 Parameters and overloads

The default behavior is the short method. Explicit configuration uses either
a parameter object or a `_with_params` suffix:

```text
mol.to_smiles()
mol.to_smiles_with_params(&params)
mol.with_3d_conformer_with_params(&params)
```

Do not create language-specific overload families with different defaults.
Bindings should expose the same explicit parameter fields and default values.

## 5. Type Classification

Public types must be classified in documentation and in the API manifest.

### 5.1 Canonical values

These are user-visible data values with value semantics:

```text
Molecule
QueryGraph
Fingerprint
DescriptorSet
Conformer
```

Detached canonical values expose accessors and local validation without
depending on parsers, matchers or runtime machinery. The live `Molecule` is
owned by the public runtime and exposes domain behavior through thin methods;
it is not subject to the detached model's dependency restriction.

### 5.2 Inputs

Parameters and options are immutable configuration values:

```text
SmilesParseParams
SmilesWriteParams
EmbedParams
MmffParams
SubstructMatchParams
```

They must have explicit defaults where the source API defines defaults. A
parameter type must not contain a live `Molecule` or a runtime capability.

### 5.3 Results and reports

Results represent calculated values, assignments, or validation reports:

```text
TopologyMapping
MatchResult
McsResult
ConformerGenerationResult
Fingerprint
```

They must not secretly contain a live `Molecule` unless the operation is
explicitly a molecule-producing public operation. Query results must remain
query results; MCS must not be forced through a concrete molecule conversion.

### 5.4 Errors

Errors are structured and stable at the public boundary. They should expose a
domain, kind, and useful detail, while allowing bindings to map them to native
exceptions or thrown errors.

Unsupported behavior must use a structured unsupported category. It must not
return an empty molecule, an empty result, a plausible placeholder value, or
panic.

### 5.5 Runtime-only values

`OpParts`, access markers, cache state, operation traces, and commit handles are
runtime implementation details. They remain private or `pub(crate)` and are
never part of the Rust facade's public API, Python classes, or WASM exports.

## 6. Cross-Language Contract

Each public API is one logical entry with three projections:

| Logical name | Rust | Python | JavaScript |
|---|---|---|---|
| `from_smiles` | `Molecule::from_smiles` | `Molecule.from_smiles` | `Molecule.fromSmiles` |
| `molecular_weight` | `mol.molecular_weight()` | `mol.molecular_weight()` | `mol.molecularWeight()` |
| `to_smiles` | `mol.to_smiles()` | `mol.to_smiles()` | `mol.toSmiles()` |
| `with_hydrogens` | `mol.with_hydrogens()` | `mol.with_hydrogens()` | `mol.withHydrogens()` |
| `tpsa` | `mol.tpsa()` | `mol.tpsa()` | `mol.tpsa()` |

The adapters may differ in:

```text
Result<T, E>       -> exception / thrown error
Option<T>          -> None / null
Vec<T>             -> list / array
&str               -> string
```

They must not differ in operation semantics, default options, result field
meaning, or behavior commitment. A declared projection specifies how a binding
must behave when implemented; it does not claim that the binding already exists.

JavaScript bindings should use `camelCase` only at the language boundary. The
Rust logical name remains the stable identifier used in manifests and tests.

## 7. API Manifest and Feature Selection

Every public method must have one manifest entry containing at least:

```text
logical name
receiver category
Rust signature
Python projection
JavaScript projection
feature capability
input and output types
error category
value-style or in-place behavior
function status
```

Fine-grained Cargo features gate capabilities, not unrelated implementation
dependencies. A feature may make its own `Molecule` method available, but it
must not expose another domain's public methods merely because the
implementation reuses an internal crate.

The user-facing `full` bundle is an explicit composition of fine-grained
capabilities and is enabled by default.

The operation registry remains the source of truth for topology mutation and
contract metadata. The public API manifest is the source of truth for naming,
receiver classification, and language projections. The two registries must
refer to the same logical operation rather than define duplicate behavior.

### Function status

Each function has one behavior status, shared by its operation metadata and
language projections rather than independently assigned in each registry:

| Status | Meaning |
|---|---|
| `Parity` | Follows the pinned upstream behavior, including options, boundary cases and errors. |
| `ParityWithDifferences` | Follows the pinned upstream except for explicitly approved, documented differences. |
| `Native` | Implements project-defined behavior with no upstream equivalence claim, such as ConfSeq. |
| `Experimental` | Offers an actual callable implementation whose behavior or API is not yet a settled commitment; limitations must be documented. |

`Parity` and `ParityWithDifferences` identify their reference library, such as
RDKit or Gemmi. `ParityWithDifferences` carries one explanation string stating
the affected conditions, the different behavior and its deliberate rationale.
It needs no separate difference ID or three-part explanation schema. Unlisted
behavior remains subject to the upstream contract.

Parity is the default development requirement for upstream ports, but the
default registry label is `Experimental` until explicitly changed. These are
different decisions: what to implement, and what commitment to publish.
Labels are maintained manually; test execution does not modify them or grant
runtime permissions. An experimental label does not relax signature checks,
operation capabilities, state validation or explicit error handling.

The registry contains real public functions and their associated public types.
Every entry must resolve to its declared Rust item; callable signatures must
be checked by the compiler. There is no separate exposure classification and
no registry placeholder for an unimplemented interface. Planned APIs belong
in the implementation plan. Data preservation is part of a function's behavior,
not another support level; structured unsupported errors remain errors, not
function status labels.

## 8. API Development Workflow

Implement each behavior once in its owning crate and expose it through the
canonical public entry point. Compatibility aliases require an explicit
decision; they must not become separate implementations.

New code must not add any of these public patterns:

```text
cosmolkit_core::Molecule
free functions taking &Molecule for ordinary molecule queries
calc_* public descriptor names
mol_to_* or mol_from_* public facade names
public OpParts or capability objects
binding-specific chemistry implementations
```

When adding or revising an API:

1. Define the logical name, signature, type classification and behavior contract
   in the implementation plan before implementing the API.
2. Declare any required operation contract before implementing its body.
3. Use the canonical public name, even when the upstream function has a
   different name. Do not automatically reproduce upstream aliases.
4. Add the `cosmolkit::Molecule` method or justified domain function together
   with its canonical binding-registry entry and compiler-checked signature.
5. Make the method extract authorized model blocks and call the algorithm
   crate.
6. Define consistent language projections. Implement bindings when they are
   within the task's scope; do not claim delivery from a projected name alone.
7. Update public documentation and examples to match the delivered API.

## 9. Review Checklist

Before adding or revising a public API, verify:

- Does the operation have a natural `Molecule` receiver?
- If yes, is it a `Molecule` method rather than a new free function?
- Is the name free of `calc_`, `mol_to_`, `mol_from_`, and redundant `get_`?
- Is value-style versus in-place behavior explicit?
- Are parameters, results, and errors separate public types?
- Does the implementation receive detached model blocks rather than `Molecule`?
- Is the operation registered with the runtime contract when it mutates state?
- Does the registry entry resolve to a real Rust item with the declared signature?
- Does the function have one behavior status and consistent language projections?
- Are feature gates scoped to the owning capability?
- Does unsupported behavior fail with a structured error?
- Is there exactly one implementation and one authoritative `Molecule`?

Examples illustrate API design, not an inventory of implemented functions.
The registry describes the actual public Rust surface. Behavior declarations
and validation results are distinct; neither substitutes for the other.

## 10. Registry Correction Guide

The four-state design above is approved but has not yet been implemented in
the registry code. Apply these corrections without changing chemical algorithms
or operation authority:

- Replace the old support/parity combinations with the single function status;
  remove `exposure` and require actual-item/signature checks for every entry.
- Keep unimplemented interfaces outside the registry. The four placeholder
  entries for `MolBlockReadParams`, `MolBlockError`, `Molecule::from_molblock`
  and `Molecule::from_molblock_with_params` have been removed.
- For this correction, mark only `fuzzy_and` and `fuzzy_or` as `Parity` against
  RDKit, on both `SparseCountFingerprint` and `SparseCountFingerprint32`.
  Initialize all other function entries as `Experimental`; do not mechanically
  preserve earlier parity claims.
- Keep status declarations independent of test execution. This correction does
  not introduce a test-to-registry promotion mechanism.
