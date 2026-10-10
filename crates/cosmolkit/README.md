# COSMolKit

`cosmolkit` is the public Rust API and runtime owner for COSMolKit 0.5.0, a Rust-native cheminformatics and structural biology toolkit. It is the sole supported Rust entrypoint for `Molecule`, operation contracts, and domain APIs. Dedicated workspace crates own detached model values and source-backed algorithms behind this facade.

The current chemistry reference is RDKit **2026.03.6** (Python distribution
`2026.3.6`), revision `0e0d85f4ca34aeae15dfc0f7cf5503bdb0a8e985`.

## Feature dependency tree

Public feature relationships.
Each edge means that enabling the parent also exposes the child's public APIs.
Default: `full`.

```text
full
├── core
├── search
│   └── core
├── reaction
│   └── search
│       └── core
├── depict
│   └── core
├── fingerprints
│   └── core
├── descriptors
│   └── core
├── conformer
│   └── core
├── stereoisomers
│   └── core
├── tautomer
│   └── search
│       └── core
├── inchi
│   └── core
└── bio
```

Top-level nodes are public feature bundles. `core` includes ordinary molecular
IO, serialization and batch processing; `serialization` and `batch` are included
capabilities, not separate user selections. There is no public `io` feature.
Batch operations for optional domains still require those domains to be enabled.

`conformer` includes 3D conformer generation, coordinate alignment, and UFF/MMFF
parameterization, energy, gradients, optimization and persistent force-field
objects. There is no separate public `forcefields` selection;
the internal conformer, force-field and alignment crates retain their
existing ownership boundaries.

The WASM distribution excludes the binary serialization path: its APIs and
implementation are not compiled, and binary-only dependencies are not included.
Other `core` capabilities, including batch processing, remain available in WASM.

Internal reuse of search or depiction helpers
does not expose those domains. IO functionality requiring query/search support
needs `search`; requesting coordinate generation needs `depict`. Dedicated APIs
are feature-gated; shared entrypoints return an explicit capability error when
the selected input or option requires a disabled domain.

For a concise Rust-native cheminformatics overview, see <https://tools.cosmol.org/rust-cheminformatics>.

## Cargo features

Default features enable `full`. Most users can keep the default or select
plain-name bundles such as `core`, `bio`, or `fingerprints`:

<!-- rust-install-version:start -->
```toml
cosmolkit = { version = "0.5.0", default-features = false, features = ["core", "bio"] }
```
<!-- rust-install-version:end -->

| Bundle | Area |
|---|---|
| `core` | Molecular text/file IO, native binary archives, batch processing, SMILES, valence, hydrogens, aromaticity, kekulization, sanitization, rings, basic stereo, matrices and transforms |
| `bio` | Structural biology readers, values, selection and operations |
| `descriptors` | Molecular descriptors |
| `tautomer` | Tautomer capability |
| `conformer` | 3D conformers, ConfSeq, alignment, UFF/MMFF energy, gradients and optimization |
| `fingerprints` | Fingerprints and molecular hashing |
| `search` | SMARTS and substructure search |
| `reaction` | SMIRKS, reaction templates and execution, including search and SMILES |
| `depict` | 2D layout and depiction |
| `inchi` | InChI and InChIKey conversion |
| `stereoisomers` | Stereoisomer enumeration |
| `full` | All features above |

`core` includes molecular text/file parsing and writing; there is no separate
`io` bundle. It does **not** include descriptors, search, depiction, tautomers or
stereoisomer enumeration. Select `descriptors`, `search`, `depict`, `tautomer`, or
`stereoisomers` explicitly when needed. Basic stereo assignment remains in
`core`; enumeration is a separate capability. Bundle names select the
documented public APIs; internal helper reuse does not expose another domain.

With defaults disabled, `features = ["bio", "core"]` does not pull in
`cosmolkit-descriptors` or `cosmolkit-tautomer` through these selections.
BIO alone activates only the structural-biology branch of IO.
Query IO requires `search`; automatic coordinate generation requires `depict`.
Ordinary molecular IO does not expose either domain. Another dependency's additive features
can still enable these packages; inspect the resolved build graph, not just
package entries in `Cargo.lock`.

For example, to enable reaction processing and its prerequisites:

```sh
cargo add cosmolkit --no-default-features --features reaction
```

With defaults disabled, neither `full` nor `core` is implicit. Features can be
combined; adding features without disabling
defaults keeps `full` enabled. Cargo features are additive across dependencies.

Each feature enables its APIs and their prerequisites. Sharing an implementation
crate does not expose unrelated functionality or promise per-function compilation.
Feature selection does not change an enabled operation's behavior or status.
`op-contracts-strict` separately enables runtime and operation-contract checks.

For exact bundle membership and registry rules, see
[Cargo feature selection](https://github.com/cosmol-studio/COSMolKit/blob/main/dev/public_api_design.md#cargo-feature-selection).

## Documentation

- [Rust API documentation](https://docs.rs/cosmolkit/latest/cosmolkit/) — types, methods, and module reference.
- [Python documentation](https://kit.cosmol.org/) — installation, API reference, and Python workflow examples.
- [Python package on PyPI](https://pypi.org/project/cosmolkit/) — package releases and installation downloads.
- [COSMolKit Web Tools](https://tools.cosmol.org/tools) — browser-based SMILES-to-SVG depiction, molecular format conversion, 3D conformer generation, InChI/InChIKey conversion, molecular property calculation, and SMILES canonicalization.
- [Validation scope and evidence](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md) — reference versions, test corpora, and documented parity boundaries.

## Rust Crates

The COSMolKit packages below are published on crates.io. Start with
`cosmolkit` for the public molecule API; the other crates separate shared values,
algorithms, and code generation. The descriptions identify each package's area;
publication does not imply that every capability in that area is implemented.
Consult the package API and support status for the selected version.

| Crate | Area |
|---|---|
| [cosmolkit](https://crates.io/crates/cosmolkit) | Public Molecule API, operation contracts, and runtime |
| [cosmolkit-types](https://crates.io/crates/cosmolkit-types) | Element, bond, and stereochemistry vocabulary |
| [cosmolkit-model](https://crates.io/crates/cosmolkit-model) | Detached atoms, bonds, topology, coordinates, properties, and query values |
| [cosmolkit-core](https://crates.io/crates/cosmolkit-core) | Foundational chemistry, ring perception, valence, and shared graph algorithms |
| [cosmolkit-macros](https://crates.io/crates/cosmolkit-macros) | Operation and binding-contract code generation |
| [cosmolkit-ringdecomposer](https://crates.io/crates/cosmolkit-ringdecomposer) | Graph cycle decomposition and Unique Ring Families |
| [cosmolkit-cx](https://crates.io/crates/cosmolkit-cx) | CX extension syntax and records |
| [cosmolkit-smiles](https://crates.io/crates/cosmolkit-smiles) | SMILES parsing, writing, canonical ranking, and CXSMILES |
| [cosmolkit-search](https://crates.io/crates/cosmolkit-search) | SMARTS queries and substructure matching |
| [cosmolkit-io](https://crates.io/crates/cosmolkit-io) | Molecular file formats and detached structure IO |
| [cosmolkit-inchi](https://crates.io/crates/cosmolkit-inchi) | InChI and InChIKey generation and conversion |
| [cosmolkit-descriptors](https://crates.io/crates/cosmolkit-descriptors) | Molecular properties and descriptors |
| [cosmolkit-fingerprints](https://crates.io/crates/cosmolkit-fingerprints) | Molecular fingerprints and similarity primitives |
| [cosmolkit-stereo](https://crates.io/crates/cosmolkit-stereo) | Stereochemistry and stereoisomer operations |
| [cosmolkit-tautomer](https://crates.io/crates/cosmolkit-tautomer) | Tautomer transformations and enumeration |
| [cosmolkit-conformer](https://crates.io/crates/cosmolkit-conformer) | 3D conformer generation and selection |
| [cosmolkit-forcefields](https://crates.io/crates/cosmolkit-forcefields) | Molecular force fields, energy, and optimization |
| [cosmolkit-alignment](https://crates.io/crates/cosmolkit-alignment) | Coordinate alignment and RMSD |
| [cosmolkit-depict](https://crates.io/crates/cosmolkit-depict) | 2D molecular layout and depiction |
| [cosmolkit-batch](https://crates.io/crates/cosmolkit-batch) | Detached batch processing and ordered results |
| [cosmolkit-bio](https://crates.io/crates/cosmolkit-bio) | Structural-biology values and hierarchy operations |

## Validation Status

**No differences have been observed on the known corpus of several hundred
thousand molecules; million-scale 0.5.0 parity validation is pending.** See
[VALIDATION.md](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md)
for the current status and the separate historical 0.3.0 evidence.

Validation compares source-defined discrete results exactly and numerical
outputs under each suite's declared comparison rule. Source reproduction,
local regressions and corpus comparisons are separate evidence boundaries.
The standard preparation and comparison workflow is documented in
[`parity-tests_fixed/README.md`](../../parity-tests_fixed/README.md).
Historical ChEMBL counts and pass results remain in
[`VALIDATION.md`](../../VALIDATION.md); they are not 0.5.0 results.

## Installation

```toml
cargo add cosmolkit
```

## Quick Start

```rust
use cosmolkit::{Molecule, SmilesWriteParams};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mol = Molecule::from_smiles("CCO")?;
    let mol = mol.with_2d_coordinates()?;

    let smiles = mol.to_smiles_with_params(&SmilesWriteParams::default())?;
    let svg = mol.to_svg(300, 300)?;

    println!("{smiles:?}");
    println!("{}", svg.len());
    Ok(())
}
```

## Molecule Operations

Normal `Molecule` operations return new values and leave the receiver
unchanged:

```rust
let mol = Molecule::from_smiles("CCO")?;
let with_h = mol.with_hydrogens()?;
assert_ne!(mol.num_atoms(), with_h.num_atoms());
```

In-place operations are explicit and always end with `_`:

```rust
let mut mol = Molecule::from_smiles("CCO")?;
mol.add_hydrogens_()?;
mol.sanitize_()?;
```

The trailing underscore is reserved for in-place mutation on public `Molecule`
methods; it has no other meaning. In-place operations prioritize avoiding the
operation-system working-copy clone when molecule blocks are uniquely owned. If
an in-place operation returns an error, the receiver is not guaranteed to equal
its pre-call value and may retain partial changes, while its internal storage
remains complete. Use the non-mutating operation when failure-preserving value
semantics are required.

Stable molecule operations include assigning atom chiral tags from a selected
3D conformer through `with_chiral_tags_from_structure()` and its explicit
in-place counterpart `assign_chiral_tags_from_structure_()`.

## Protein Structures

```rust
use cosmolkit::Protein;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string("1crn.pdb")?;
    let protein = Protein::from_pdb(&text)?;
    let summary = protein.selection_summary();

    println!("chains: {}", summary.chains);
    println!("residues: {}", summary.residues);
    println!("atoms: {}", summary.atoms);
    Ok(())
}
```

## Batch Workflows

```rust
use cosmolkit::MoleculeBatch;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = vec![
        "CCO".to_string(),
        "c1ccccc1".to_string(),
        "CC(=O)O".to_string(),
    ];

    let batch = MoleculeBatch::from_smiles_list(&smiles)?
        .with_parallel_jobs(Some(8))?
        .with_2d_coordinates()?;

    let out = batch.to_smiles_list()?;
    println!("{out:?}");
    Ok(())
}
```

## InChI

Enable `inchi` alongside `core` (or use `full`) for molecule conversion and
direct InChIKey generation:

```rust
use cosmolkit::{Molecule, inchi_to_key};

let molecule = Molecule::from_smiles("C")?;
let identifier = molecule.to_inchi()?;
assert_eq!(identifier, "InChI=1S/CH4/h1H4");

let key = inchi_to_key(&identifier)?;
assert_eq!(key, "VNWKTOKETHGBQD-UHFFFAOYSA-N");

let parsed = Molecule::from_inchi(&identifier)?;
assert_eq!(parsed.to_inchi()?, identifier);
```

`InchiReadParams` controls sanitization and hydrogen removal (both default to
true). `InchiWriteParams` carries the engine option string. The corresponding
`*_with_params` methods accept these values. Queries leave the molecule
unchanged; failures return typed `InchiError` values. The implementation uses
the official InChI v1.07.5 engine and RDKit 2026.03.6 adapter source. IXA, AuxInfo
reconstruction, INCHIGEN, version-query, and extended-polymer entry points are
not exposed by this facade.

## Molecular Descriptors

Enable `descriptors` (or use `full`). Public molecule queries delegate to the
detached `cosmolkit-descriptors` algorithms:

```rust
use cosmolkit::Molecule;

let molecule = Molecule::from_smiles("c1ccccc1O")?;
assert_eq!(molecule.molecular_formula()?, "C6H6O");
assert!(molecule.molecular_weight()? > 94.0);
assert_eq!(molecule.num_aromatic_rings()?, 1);
assert!(molecule.chi_0()? > 0.0);
assert_eq!(molecule.mqns(false)?.len(), 42);
```

The descriptor surface includes molecular properties, connectivity and shape
indices, Lipinski and ring/stereo counts, MQN, Labute ASA, and SlogP/SMR VSA.
Queries use their declared prepared-state requirements and return typed
`DescriptorReadError` values; they do not silently install missing caches.

```rust
let mol = Molecule::from_smiles("CCO")?;
assert_eq!(mol.num_heavy_atoms()?, 3);
assert_eq!(mol.total_atom_count()?, 9);
assert_eq!(mol.lipinski_hba()?, 1);
assert_eq!(mol.lipinski_hbd()?, 1);
assert_eq!(mol.fraction_csp3()?.to_bits(), 1.0_f64.to_bits());
```

`num_atoms()` counts explicit atom rows; `total_atom_count()` also includes
implicit and explicit-property hydrogens. `lipinski_hba()` is the N/O count,
and `lipinski_hbd()` is the donor-hydrogen sum on N/O, not the general
acceptor/donor-atom SMARTS counts.

## Fingerprints

Enable `fingerprints` (or use `full`). The public `fingerprint_*` family
covers Morgan, MACCS, AtomPair, Topological Torsion, RDKit topological, Pattern,
Layered and Avalon fingerprints. Topological Torsion is the ordered atom-path
family, distinct from the path/subgraph Topological fingerprint.

```rust
use cosmolkit::{FingerprintAdditionalOutput, Molecule, TopologicalTorsionFingerprintParams};

let molecule = Molecule::from_smiles("c1ccccc1O")?;
let morgan = molecule.fingerprint_morgan()?;
let avalon = molecule.fingerprint_avalon()?;
let layered = molecule.fingerprint_layered()?;
let params = TopologicalTorsionFingerprintParams::default();
let mut output = FingerprintAdditionalOutput::default();
output.allocate_atom_to_bits();
output.allocate_bit_paths();
let torsion = molecule.fingerprint_topological_torsion_with_params(
    &params, Some(&mut output),
)?;
assert_eq!(morgan.n_bits(), 2048);
assert_eq!(avalon.n_bits(), 512);
assert_eq!(layered.n_bits(), 2048);
assert_eq!(torsion.n_bits(), 2048);
```

Use `*_with_params` for explicit settings and the family's supported output
collector for atom/bit/path provenance. Batch APIs preserve input order and
record positions for failures. Current corpus comparisons use RDKit 2026.03.6;
historical full-corpus totals remain in [VALIDATION.md](../../VALIDATION.md).

Layered retains the upstream algorithm's `0.7.0` metadata; this is not the
COSMolKit package version. Unspecified roots enumerate the whole molecule,
whereas an explicitly empty root list enumerates no paths. Invalid parameters
return typed errors. Native reference crashes are reported separately, not
counted as matches.

## Conformer Generation And Force Field Applications

Enable `conformer` (or use `full`) for embedding, alignment and UFF/MMFF.
Embedding returns new molecule values and supports fixed seeds, RMS pruning
and sequential seed expansion.

```rust
use cosmolkit::{EmbedParams, Molecule};

let molecule = Molecule::from_smiles("CC(=O)NC")?.with_hydrogens()?;
let mut params = EmbedParams::etkdg_v3();
params.random_seed = 123;
params.num_threads = 1;
let embedded = molecule.with_3d_conformer_with_params(&params)?;
let multiple = molecule.with_3d_conformers_with_params(5, &params)?;
assert_eq!(embedded.conformers_3d().len(), 1);
assert!(!multiple.conformers_3d().is_empty());
```

Force-field calls require existing coordinates and do not add hydrogens or
embed automatically. Value-style optimization leaves the input unchanged:

```rust
use cosmolkit::{MmffOptimizationParams, UffOptimizationParams};

if embedded.uff_has_all_molecule_params()? {
    let result = embedded.with_uff_optimized_with_params(&UffOptimizationParams::default())?;
    println!("UFF energy: {}", result.energy);
}
if embedded.mmff_has_all_molecule_params()? {
    let result = embedded.with_mmff_optimized_with_params(&MmffOptimizationParams::default())?;
    println!("MMFF needs_more: {}", result.needs_more());
}
```

For repeated interactive evaluation, use an owned persistent handle:

```rust
use cosmolkit::ForceFieldMinimizeParams;

let mut field = embedded.uff_force_field()?;
let initial_energy = field.energy()?;
let gradient = field.gradient()?;
let outcome = field.minimize_with_params_(&ForceFieldMinimizeParams::new(20, 1e-4, 1e-6))?;
println!("{initial_energy} {gradient:?} {outcome:?}");
```

`mmff_force_field()` provides the corresponding MMFF handle.
`set_positions_()` and `set_fixed_atoms_()` update the detached evaluator;
they do not implicitly write coordinates back to the source molecule.

## Molecular Alignment And RMSD

Ordinary molecular alignment uses the source-backed RDKit MolAlign boundary.
Transform queries, best RMSD, coordinate-frame RMSD, and all-conformer pair
measurements are read-only. Coordinate changes are exposed only through a
value-style method or an explicit trailing-underscore method.

```rust
use cosmolkit::{AlignmentAtomMap, AlignmentParameters, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference = Molecule::from_smiles("CCC")?.with_only_3d_conformer(
        vec![vec![0.0, 0.0, 0.0], vec![1.0, 0.0, 0.0], vec![0.0, 2.0, 0.0]],
    )?;
    let probe = Molecule::from_smiles("CCC")?.with_only_3d_conformer(
        vec![vec![3.0, -2.0, 1.0], vec![4.0, -2.0, 1.0], vec![3.0, 0.0, 1.0]],
    )?;
    let params = AlignmentParameters {
        atom_map: Some(
            (0..3)
                .map(|index| AlignmentAtomMap {
                    probe_atom: index,
                    reference_atom: index,
                })
                .collect(),
        ),
        ..Default::default()
    };

    let measured = probe.alignment_transform_to_with_params(&reference, &params)?;
    let (aligned, applied) = probe.with_alignment_to_with_params(&reference, &params)?;
    assert_eq!(probe.conformers_3d()[0].coordinates()[0], [3.0, -2.0, 1.0]);
    assert!(measured.rmsd < 1.0e-8 && applied.rmsd < 1.0e-8);
    assert_eq!(aligned.conformers_3d()[0].coordinates()[0], [0.0, 0.0, 0.0]);
    Ok(())
}
```

Weighted and reflected alignment, automatic or explicit atom maps, conformer
IDs, iteration limits, best-map selection, and conformer-set alignment use
typed parameter objects. O3A and MMFF/Crippen scoring are separate capabilities
and are not implied by this ordinary MolAlign API.

## Examples

```bash
cargo run -p cosmolkit --example smiles_write_options
cargo run -p cosmolkit --example draw_svg
cargo run -p cosmolkit --example draw_png
cargo run -p cosmolkit --example sdf_to_smiles
cargo run -p cosmolkit --example protein_from_pdb
cargo run -p cosmolkit --example read_xyz
cargo run -p cosmolkit --example molalign_rmsd
cargo run -p cosmolkit --example conformer_generation
cargo run -p cosmolkit --example forcefield_optimization
```

## Contributor Validation

Detached core regressions and facade operation contracts are separate checks:

```bash
cargo check -p cosmolkit-core
cargo test -p cosmolkit-core --release
cargo check -p cosmolkit --features op-contracts-strict
cargo test -p cosmolkit --profile dev-test --features op-contracts-strict
cargo check -p cosmolkit-py
cargo fmt --all
```

Use debug-profile test filters for small local iterations. Use the optimized
`dev-test` profile for daily suites and `release` for distribution builds.
Enable `cosmolkit/op-contracts-strict` explicitly for operation-contract checks;
optimization profiles do not enable strict checks themselves. Core has no
runtime-check features. See the repository development manual for the full
pre-commit checklist.

Python binding validation:

```bash
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
.venv/bin/pytest
```

COSMolKit 0.5.0 uses the split-crate architecture: `cosmolkit` owns the public
`Molecule` API and operation runtime; model and domain crates own detached
values and algorithms. Users import the public chemistry API from `cosmolkit`,
not implementation crates. See the repository architecture and public API
design documents for ownership, contracts and feature selection.
