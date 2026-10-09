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
| Query construction | facade function and thin value factory | `parse_smarts`; Python `QueryGraph.from_smarts` |
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
parse_smarts(text) -> QueryGraph
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

Python business functions are exposed directly under `cosmolkit`, not under
domain submodules such as `cosmolkit.search` or `cosmolkit.confseq`. Rust domain
modules remain valid implementation/semantic groupings; their layout does not
introduce an extra Python namespace. Use a domain-explicit function name when
flattening would make a generic verb ambiguous, for example `decode_confseq`
and `decode_confseq_batch`, not top-level `decode` and `decode_batch`.

The user-approved SMARTS construction surface deliberately provides both:

```python
query = cosmolkit.parse_smarts(text)
query = cosmolkit.QueryGraph.from_smarts(text)
```

Explicit configuration uses `parse_smarts_with_params(text, params)` and
`QueryGraph.from_smarts_with_params(text, params)`. Both forms return the same
canonical query value and share the same parser, defaults and typed errors.
The bound class factory is a thin projection of a registered facade callable;
it does not add a parser dependency or inherent parsing implementation to the
detached model type. Do not retain the previous `cosmolkit.search.*` layer as
a compatibility namespace. Other search functions (`compile_query`,
`write_smarts`, `write_cx_smarts`) are likewise top-level Python functions.

The top-level public `BioStructure` and `Protein` use associated constructors
`from_pdb`, `from_pdb_with_params` and `from_mmcif`. Their detached data stays
in BIO and their parsers stay in IO; do not replace these constructors with
model-prefixed module functions or add inherent IO methods to a foreign type.
Both objects' in-place operations use the trailing underscore and their
value-style counterparts leave the input unchanged, through the single
lightweight BIO operation declaration described in [BIO architecture](./bio_architecture.md).

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

Use `write_*` when the method writes a file or directory, rather than returning
the serialized representation. In particular, `Molecule::to_sdf()` and
`SdfRecord::to_sdf()` return text; `MoleculeBatch::write_sdf(path)` and
`MoleculeBatch::write_sdf_files(directory)` write files and return an export
report. Their explicit Rust configuration forms use `_with_params`; Python
and JavaScript retain the same `write_*` / `write*` semantic names.

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

### 4.6 Domain-explicit public type names

A public type name must identify its domain and role at its actual exposure
location. Top-level facade exports must remain understandable without a nearby
algorithm call, source-library knowledge, or an implementation-module name.
Generic names such as `AdditionalOutput`, `Options`, or `Result` are not
sufficient for domain-specific top-level types. Use a domain qualifier when the
containing public namespace does not already make the meaning unambiguous;
do not add redundant prefixes to already self-explanatory vocabulary types.

Counterexample:

```python
output = ck.AdditionalOutput()  # Additional output for which functionality?
```

Canonical name:

```python
output = ck.FingerprintAdditionalOutput()
fp = mol.fingerprint_morgan_with_params(params, output)
```

`FingerprintAdditionalOutput` names the optional fingerprint metadata collector
shared by multiple fingerprint algorithms, not a Morgan-only result and not a
container that owns the returned fingerprint. Keep this type name consistent
across Rust, Python, and JavaScript. RDKit's original `AdditionalOutput` spelling
remains verbatim in pinned-source anchors; upstream names do not override the
public naming rule. Register and update the canonical type, signatures, bindings,
documentation, and tests together. Do not retain the ambiguous name as a
compatibility alias without explicit authorization.

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

Parameters and options are configuration values. Their Rust API is unchanged
by the language-projection rules below:

```text
SmilesParseParams
SmilesWriteParams
EmbedParams
MmffParams
SubstructMatchParams
```

They must have explicit defaults where the source API defines defaults. A
parameter type must not contain a live `Molecule` or a runtime capability.

#### Python and JavaScript: writable configuration properties

Python and JavaScript parameter/option objects must support both constructor
configuration and assignment after default construction. Every public input
field must have a getter and setter, not only selected fields or types:

```python
params = ck.SubstructMatchParams()
params.use_chirality = True
params.max_matches = 100
# Equivalent configuration:
params = ck.SubstructMatchParams(use_chirality=True, max_matches=100)
```

```javascript
const params = new ck.SubstructMatchParams();
params.useChirality = true;
params.maxMatches = 100;
```

This requirement applies only to Python and JavaScript configuration objects;
it does not make molecule storage, calculated results, reports, or errors
writable, and does not change Rust APIs or molecule operation contracts.
Computed properties are not input fields and may remain read-only.

Constructors and setters must use the same field types and validation rules.
An invalid assignment raises the language's corresponding error and leaves
the previous field value unchanged. Changing a configuration object must not
retroactively change an earlier operation or another independent object.
Register writable properties consistently with the binding contract, generate
matching `.pyi` and TypeScript declarations, and test constructor/assignment
equivalence through an actual operation. A read-only binding for a configurable
input field is a contract defect, not an alternative configuration style.

Parameterized operations must also expose one idiomatic method name supporting
default configuration, a parameter instance, or convenient field configuration:

```python
mol.substruct_matches(query)
mol.substruct_matches(query, params)
mol.substruct_matches(query, use_chirality=True, uniquify=False, max_matches=100)
```

```javascript
mol.substructMatches(query);
mol.substructMatches(query, params);
mol.substructMatches(query, {
    useChirality: true,
    uniquify: false,
    maxMatches: 100,
});
```

Python configuration keywords are keyword-only and mutually exclusive with a
parameter instance; supplying both raises `TypeError`, without implicit
overrides. JavaScript accepts either a parameter instance or a plain options
object as the configuration argument; it does not have Python-style keyword
arguments. Omitted configuration fields use the same registered defaults as
the parameter constructor. Reject unknown fields and invalid values rather
than silently ignoring them.

Bindings normalize these forms to the same Rust parameter type and canonical
operation; they must not duplicate algorithms or change Rust signatures.
Python `.pyi` files express the call forms using `@overload`; TypeScript
declarations expose the parameter-instance/options-object forms. TypeScript
implementations may define default values, but `.d.ts` declarations describe
optional parameters/properties without default-value initializers. Document
the defaults and test all call forms for equivalent behavior. A `_with_params`
suffix must not be required to select explicit configuration in Python or
JavaScript.

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

### Enforced projection boundary

Binding delivery is gated by the compiled canonical registry, not a manually
maintained list or an author's checklist. The stub generator and WASM package
validator must reject a projection that violates any of these requirements:

1. A `parameter` record has a registered configuration constructor schema.
   Rust public-field records may use `Default` and struct construction; a
   redundant Rust `new()` is not required. Language-specific constructor fields,
   types and defaults are declared on that type in the canonical registry and
   checked against compiled binding metadata and actual descriptors. Its constructor
   inputs define the public configurable fields; each must be readable and
   writable in Python and JavaScript declarations and actual descriptors.
   Constructor, getter and setter types must agree. Computed output fields are
   not configurable inputs. Enums, bit masks and opaque selectors use the
   explicit `parameter_selector` role, not a name-based exemption.
   A Python constructor is registered on its type, for example
   `python_configuration: [{ name: max_matches, python_type: "builtins.int", default: "1000" }]`.
   Defaults use Python expressions; `default: required` denotes a required
   argument. This describes the language constructor, not a fictitious Rust
   callable. The checker compares this schema with compiled constructor
   declarations and exercises real properties, including nested configuration
   preservation and failure atomicity.
2. Canonical default/configured operation pairs expose the same language method
   with parameter-object and field-based call forms. Field configuration uses
   the constructor's field types and defaults. Python declares explicit
   overloads and keyword-only configuration; TypeScript must accept the actual
   parameter instance and the options object. Bindings reject mixed/unknown
   configuration at runtime rather than discarding it.
3. Every enabled registry type, callable and declared property exists in the
   final language surface. Check the actual Python module and generated WASM
   exports as well as `.pyi`/`.d.ts`. A binding source and its generated output
   both omitting an API is a failure. Native binary archive exclusion on WASM
   follows the platform contract; missing implementations and Experimental
   status are not exemptions.
   A scalar value may explicitly declare a native Python projection such as
   `builtins.int`; this does not promise an independently exported wrapper class.
   Configuration records and object adapters cannot use that scalar projection.

`stub_gen` runs these checks before replacing the existing `.pyi`; the WASM
build runs them before exporting a distribution. Neither gate may synthesize
declarations for missing implementations. Existing violations remain failed
delivery gates until the bindings or their incorrect registration are fixed.
Author review supplements these structural constraints. Small actual-operation
regressions still prove parameter forwarding and failed-assignment atomicity;
declaration checks alone cannot prove either behavior.

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

The operation registry remains the source of truth for topology mutation and
contract metadata. The public API manifest is the source of truth for naming,
receiver classification, and language projections. The two registries must
refer to the same logical operation rather than define duplicate behavior.

### Cargo feature selection

Public documentation and examples use plain feature names. Ordinary users must
not need to understand `cap-*` gates to select functionality. The crate README's
dependency tree shows user features and included functionality, not internal
selectors. `cap-*` names remain implementation-level capability gates used by
cfg conditions and registries; their details belong in development documentation.

Domain names without a prefix select user bundles; `cap-*` names select individual capabilities.
Bundles compose capabilities and their public prerequisites. Both forms can be
combined, and Cargo adds their selections together. `full` is the default.

| Bundle | Direct feature membership (prerequisites are transitive) |
|---|---|
| `core` | `cap-smiles`, `cap-hydrogens`, `cap-valence`, `cap-radicals`, `cap-rings`, `cap-matrices`, `cap-transforms`, `cap-stereo`, `cap-kekulize`, `cap-aromaticity`, `cap-sanitize`, `cap-io`, `cap-serialization`, `cap-batch` |
| `bio` | `cap-bio` |
| `descriptors` | `cap-descriptors`, `core` |
| `tautomer` | `cap-tautomer`, `search` |
| `conformer` | `cap-conformer`, `cap-confseq`, `cap-alignment`, `cap-forcefields`, `core` |
| `fingerprints` | `cap-fingerprints`, `cap-hashing`, `core` |
| `search` | `cap-search`, `core` |
| `reaction` | `cap-reaction`, `search` |
| `stereoisomers` | `cap-stereoisomers`, `core` |
| `depict` | `cap-depict`, `core` |
| `inchi` | `cap-inchi`, `core` |
| `full` | All bundles above |

Molecular format parsing and writing, whether from strings or files, belong to
the `core` bundle; there is no top-level `io` bundle. The internal IO crate
remains the unique implementation owner. Native binary archives and batch
processing are included in `core`, without separate plain-name selectors.
WASM excludes binary archive APIs, implementation and binary-only dependencies.

SMIRKS parsing, reaction templates and the public `Reaction` API use
`cap-reaction`, selected by the `reaction` bundle and included in `full`.
`reaction` includes `search`, which includes `core`: callers do not need to
select those prerequisites separately to use their public APIs.

`core` means foundational molecule chemistry and SMILES, not all inexpensive
or historical APIs. Descriptors, search, depiction, tautomer capability and
stereoisomer enumeration require explicit selection outside `core`. Basic
stereo assignment is part of `core`; `stereoisomers` selects enumeration
separately and remains included in `full`. Feature membership does not claim
that every API planned for the domain is implemented.

Most callers use the default or select groups such as `core` and `bio`.
For precise selection, disable defaults and choose individual capabilities:

```sh
cargo add cosmolkit --no-default-features --features cap-io,cap-kekulize,cap-sanitize,cap-hydrogens
```

With defaults disabled, `core` is not implicit. An empty selection retains
the live molecule/runtime, builders and foundational model values, but no
optional capability APIs. Adding `features` without disabling defaults keeps
`full` enabled. Cargo feature unification also means another dependency can
enable more capabilities; features are additive, not a deny list.

Public cfg gates and both registry feature fields use the owning `cap-*`
selector. Always-present declarations use `runtime` or `metadata` labels;
those labels are not Cargo capability selectors. Bundle membership is defined
in `crates/cosmolkit/Cargo.toml`, not duplicated in a production registry.
The sole platform-gate exception is native binary serialization:
`all(feature = "cap-serialization", not(target_arch = "wasm32"))`.
It excludes native archives under WASM `full` without introducing another
public selector. Other compound registry gates remain forbidden.

A public bundle exposes its declared prerequisites: reaction and tautomer
include search; every molecular bundle includes core. BIO remains independent.
The public conformer bundle combines generation, alignment and forcefields,
without changing their internal owners. Private search, stereo, forcefield or
depiction helpers do not by themselves expose another public domain.
Query IO requires search; automatic coordinate generation requires depict.
Dedicated APIs are gated and shared entrypoints return explicit capability
errors when the requested behavior is unavailable. Do not silently discard
query semantics or stereochemistry to avoid a dependency.
The dependency tree belongs at the top of `crates/cosmolkit/README.md` and must
match the Cargo manifest. WASM's feature table is generated from that manifest,
with only documented platform exclusions, not a second hand-maintained tree.

Sharing a foundational implementation crate does not select all algorithms in
that crate. For example, `cap-smiles` uses `cosmolkit-core` but does not select
every capability in the `core` bundle. Public prerequisites do not alter
algorithm ownership or operation authority.

IO's internal `molecule` feature gates molecular formats; `bio` independently
gates BIO formats/CID. The facade disables IO defaults: `cap-bio` selects only
the structural BIO branch, not `Molecule::from_sdf` or molecular search. The
domain IO crate defaults to all molecular format subfeatures for direct owner
builds. Its internal `search`, `depict` and native `binary` features can be
selected separately from ordinary `molecule` IO.
An isolated external consuming build selecting only `bio` and/or `core` must
not resolve descriptors or tautomer. Additional additive selections may
legitimately enable other capabilities. Check active resolved dependency edges
rather than lockfile package presence. These switches do not promise
per-function compilation inside an implementation crate.

`runtime-invariants`, `op-contracts`, and `op-contracts-strict` are separate
validation switches, not chemistry bundles.
Feature selection changes availability, not an enabled function's semantics,
operation authority or `FunctionStatus`.

### Instance receiver ownership

The binding registry declares `receiver: shared`, `receiver: mutable`, or
`receiver: owned`, independently of the result type. An omitted instance
receiver defaults to `mutable` for `in_place`, otherwise `shared`; consuming
receivers must explicitly declare `owned`. Static/module entries have no receiver.
The generated metadata exposes this as `BindingCallableContract.receiver`.

Shared receivers use `&Self`; mutable receivers use `&mut Self` and require
`in_place`. Owned receivers use `Self`, require `value_returning`, and cannot
link to the borrowed/in-place operation lifecycle. The compiler checks the
declared function signature against the actual method. No business type or
method name grants consuming permission. `into_*` conversions transfer ownership;
the trailing `_` remains reserved for in-place mutation, not consumption.

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

Declarations use `status: parity("RDKit")`,
`status: parity_with_differences("Gemmi", "approved conditions, behavior and rationale")`,
`status: native`, or `status: experimental`. Omission defaults to `Experimental`.
Reference names and approved-difference explanations must be nonempty.
For registered Molecule operations, declare the status once in `molecule_ops!`;
binding entries inherit it through the generated value/in-place method metadata
and cannot override it. Feature metadata describes capability selection only.

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
- Does each public type name identify its domain and role at its exposure location?
- Is value-style versus in-place behavior explicit?
- Are parameters, results, and errors separate public types?
- Do Python/JavaScript configuration input fields support both construction and assignment?
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

## 10. Naming and usability corrections

This section records approved, targeted corrections. Update the canonical Rust
API, registry, Python/JavaScript projections, generated declarations,
documentation and tests together. It does not authorize compatibility aliases
or changes to chemical behavior and does not imply that a listed implementation
gap has been completed.

### 10.1 Specific configuration naming corrections

The following are individual naming corrections, not a general cross-domain
spelling rule. The preference for `remove_hs` and omission of `do_` apply to
the listed cases only; do not extrapolate them to other APIs. These corrections
do not rename value/in-place molecule operations such as `without_hydrogens()`
and `remove_hydrogens_()`.

| Counterexample | Proposed positive example | Reason |
|---|---|---|
| BIO uses `remove_hs`, but SMILES/SDF/Mol2/InChI use `remove_hydrogens` | All corresponding read/conversion options use `remove_hs`; JavaScript uses `removeHs` | Users should not relearn the same option for every format. |
| A boolean is named `do_kekulize` | `kekulize=True` | `do_` adds no information. |
| SMILES uses `do_kekule` and `do_isomeric_smiles` | `kekule=True`, `isomeric_smiles=True` | Name the requested representation directly. |
| Related conformer APIs mix `confs` and `conformers` | `clear_conformers`, `with_uff_optimized_conformers`, `with_mmff_optimized_conformers` | Use the same term for the same object. |
| SDF entry points mix `coordinate_dim="auto"` and a typed `coordinate_mode` | `coordinate_mode=ck.SdfCoordinateMode.Preserve` on all corresponding entry points | Share one typed policy and its defaults. |

`kekulize` describes a preprocessing operation; `kekule` describes a requested
SMILES representation. Removing `do_` does not make those two behaviors
interchangeable. Likewise, do not rename a genuine coordinate dimension into a
coordinate policy merely because the words look similar.

### 10.2 Avoid redundant conversion names

| Counterexample | Proposed positive example |
|---|---|
| `ck.inchi_to_inchi_key(inchi_text)` | `ck.inchi_to_key(inchi_text)` |
| A generic top-level `ck.to_key(text)` | `ck.inchi_to_key(inchi_text)` |

Keep `mol.to_inchi_key()` for conversion from a molecule. The free function's
input is already InChI text; the molecule method's input is a molecule. A naming
cleanup must preserve that distinction and must not add another chemical
conversion to the text-to-key function.

### 10.3 Fingerprint-first method families

Put the shared functionality first, so typing `mol.fingerprint_` discovers the
available algorithms. Use the order **fingerprint → algorithm → output form**.
Apply it consistently to scalar, query and batch entry points.

| Counterexample | Proposed positive example |
|---|---|
| `mol.layered_fingerprint()` | `mol.fingerprint_layered()` |
| `mol.pattern_fingerprint()` | `mol.fingerprint_pattern()` |
| `mol.morgan_fingerprint()` | `mol.fingerprint_morgan()` |
| `mol.atom_pair_fingerprint()` | `mol.fingerprint_atom_pair()` |
| `mol.topological_fingerprint()` | `mol.fingerprint_topological()` |
| `mol.topological_torsion_fingerprint()` | `mol.fingerprint_topological_torsion()` |
| `mol.maccs_fingerprint()` | `mol.fingerprint_maccs()` |
| `mol.morgan_count_fingerprint()` | `mol.fingerprint_morgan_count()` |
| `mol.morgan_sparse_fingerprint()` | `mol.fingerprint_morgan_sparse()` |
| `mol.morgan_sparse_count_fingerprint()` | `mol.fingerprint_morgan_sparse_count()` |
| Scalar `layered_fingerprint`, but batch `fingerprint_layered_list` | Scalar `fingerprint_layered`, batch `fingerprint_layered_list` |
| `layered_query_fingerprint(...)` | `fingerprint_layered_query(...)` |

This is a discoverability rule for entry points, not permission to collapse
different algorithms, bit/count representations, sizes or legacy semantics.
Do not apply redundant prefixes to a generator's own `fingerprint()` method
when the receiver already identifies the functionality.

### 10.4 File writes must read as file writes

| Counterexample | Proposed positive example |
|---|---|
| `batch.to_images(directory)` creates directories and writes files | `batch.write_images(directory)` |
| `batch.to_images_with_params(directory, params)` writes files | Rust `batch.write_images_with_params(directory, &params)`; Python `batch.write_images(directory, params)` |

Keep `mol.to_svg()` for returned SVG text and `mol.write_svg(path)` for a file.
The same conversion/write distinction already applies to SDF; images are not
an exception.

### 10.5 Configuration must not create a second Python API vocabulary

Counterexample:

```python
mol.layered_fingerprint_with_output_with_params(params)
# Users must discover a separate method name just to supply configuration.
```

Proposed positive examples:

```python
mol.fingerprint_layered_with_output()
mol.fingerprint_layered_with_output(params)
mol.fingerprint_layered_with_output(fp_size=1024)

mol.substruct_match(query, use_chirality=True)
mol.has_substruct_match(query, use_chirality=True)
mol.substruct_matches(query, use_chirality=True)

batch.write_images(
    directory,
    ck.BatchImageParams(execution=ck.BatchParams(errors="keep")),
)
```

The last example preserves the distinction between image options and batch
execution options; it does not invent a second error policy. Rust may retain
explicit `_with_params` methods. Python/JavaScript expose the same operation
through default, parameter-object and keyword/options-object forms, with the
same defaults and errors. Do not make configurable single/boolean substructure
queries default-only while the multiple-match query accepts configuration.
Generated declarations must describe working calls, not desired overloads
that the runtime binding rejects.

### 10.6 Cache observation must be explicit

Counterexample: document `mol.atoms()` as merely inspecting existing cached
valence when it actually requests recalculation, or use it in a test intended
to observe unsanitized cached values.

Positive examples using the current explicit metadata API:

```python
cached = mol.atom_metadata(recalculate=False)  # Observe existing valid cache.
fresh = mol.atom_metadata(recalculate=True)   # Request recalculation.
```

Document that `atoms()` currently requests recalculated metadata. A cache-only
read must not silently recalculate, and a recalculating read must not be
advertised as cache-only. This clarification does not approve changing cache
validity, failure semantics or operation permissions.

### 10.7 Review names separately from implementation gaps

Counterexample: report a feature as missing solely because its previous public
name raises `AttributeError`.

Positive example: distinguish the old
`ck.get_topological_torsion_fingerprint_as_ids(mol)` call from the existing
`mol.topological_torsion_ids()` method, then verify the actual returned IDs.

Classify findings as a naming/test update, a missing binding, missing facade
integration, an incomplete algorithm branch, or an unresolved behavior
difference. Do not solve an implementation gap with an alias, fabricate a stub,
or change chemical expectations to match the implementation.


### 10.8 Python: enum values and string inputs

Every public Python input whose logical type is an enum must accept both an
enum member and its documented string spelling. This applies to function and
method arguments, configuration constructors, and writable configuration
properties, not only to `coordinate_mode`.

```python
params = ck.SdfReadParams(coordinate_mode=ck.SdfCoordinateMode.Preserve)
params = ck.SdfReadParams(coordinate_mode="preserve")
params.coordinate_mode = "require_3d"
params.coordinate_mode = ck.SdfCoordinateMode.Require3D

mol = ck.Molecule.from_sdf(text, coordinate_mode="preserve")
mol = ck.Molecule.from_sdf(text, coordinate_mode=ck.SdfCoordinateMode.Preserve)
```

For `SdfCoordinateMode`, the canonical string spellings are `"preserve"`,
`"require_2d"` and `"require_3d"`. Define each enum's string vocabulary from
its declared members and document it; do not guess values by fuzzy matching,
silently fall back to a default, or introduce different vocabularies for
different entry points. Equivalent enum and string inputs must normalize to
the same Rust enum and use the same validation and algorithm.

Native Python enum spellings are generated deterministically from declared
variant names using lowercase snake case (`NonStrict` → `"non_strict"`,
`Require3D` → `"require_3d"`, `V2000` → `"v2000"`). The same declaration
generates extraction and the vocabulary consumed by binding checks; do not
maintain separate per-function conversion tables.

Unknown strings raise `ValueError`; values of an unrelated type raise
`TypeError`. Failed property assignment leaves the previous value unchanged.
An optional enum accepts `None` only when the underlying contract is optional.
Getters and results retain their canonical enum type; accepting string inputs
does not convert enum outputs into strings or change Rust signatures.

Generate input annotations such as `SdfCoordinateMode | str` in `.pyi`, while
keeping getter/result annotations as `SdfCoordinateMode`. Registry-backed
binding checks and focused Python regressions must verify constructor,
assignment and direct-call equivalence, invalid-input rejection and unchanged
defaults. A stub advertising string support without an actual working binding
does not satisfy this requirement.
