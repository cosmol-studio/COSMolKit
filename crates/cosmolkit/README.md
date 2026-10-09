# COSMolKit

`cosmolkit` is the public Rust API and runtime owner for COSMolKit 0.5.0, a Rust-native cheminformatics and structural biology toolkit. It is the sole supported Rust entrypoint for `Molecule`, operation contracts, and domain APIs. Dedicated workspace crates own detached model values and source-backed algorithms behind this facade.

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
cosmolkit = { version = "0.5.0-rc.19", default-features = false, features = ["core", "bio"] }
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
`core`; enumeration is a separate capability. Bundle names select features,
not a promise that every planned API in that area is already implemented.

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

**0.5.0 validation is pending.** The results summarized below are from
**0.3.0**, not a validation pass for 0.5.0. See
[VALIDATION.md](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md)
for the historical boundary and the current pending status.

### Historical 0.3.0 evidence

COSMolKit treats parity as **source-backed semantic equivalence within explicitly documented boundaries**, not as statistical agreement of final outputs. Compatibility-critical chemistry is implemented as a line-by-line, source-backed port with explicit operation contracts and traceable correspondence to pinned upstream code. Validation corpora verify that port; they are not used to iteratively tune heuristic reimplementations until outputs happen to agree.

The comparison boundary therefore extends well beyond final strings. Covered surfaces compare exact bytes, bits, return status, complete atom and bond state, stereochemistry, derived state and invariants, **RNG state, seed handling, and random draw sequences where stochastic behavior is part of the contract**, every matrix entry, coordinates, energies, and every gradient component where applicable. Discrete results must match exactly; declared numerical tolerances reach `1e-8` for matrix entries and `1e-6` for coordinates, energies, and gradients. **99% or 99.9% agreement remains unfinished when any covered mismatch exists.**

This boundary is stress-tested against a complete ChEMBL 37 profile: 2,897,819 source records, 2,897,804 of them mutually parseable, across 34 repository-defined sharded phases against pinned RDKit `2026.03.1`. The profile performs billions of comparisons, expands parameter spaces into matrices of up to 768 branches, repeats complete matrices to expose instability, permutes operation order, and checks scalar, one-thread, multi-thread, batch, and shared-object concurrent paths.

Every discovered mismatch is traced back to the corresponding upstream logic, corrected at the source-port level, and permanently retained as a focused regression rather than hidden by corpus-specific adjustments. This discipline limits **semantic debt** by preventing convenient local fixes from accumulating into undocumented chemistry behavior.

The parity suite uses three complementary validation layers. The complete ChEMBL 37 profile provides large-scale stress coverage; the maintained 5,000-record corpus runs exhaustive parameter matrices not yet practical across the full ChEMBL profile; and the 152-record project corpus keeps focused regressions fast enough for daily testing.

See [`VALIDATION.md`](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md) for exact corpus eligibility, comparison counts, tolerances, per-feature boundaries, focused regressions, and upstream surfaces outside the current claim.

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

    println!("{smiles}");
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
    let protein = Protein::from_pdb("1crn.pdb")?;
    let summary = protein.selection_summary();

    println!("chains: {}", summary.chains);
    println!("residues: {}", summary.residues);
    println!("atoms: {}", summary.atoms);
    Ok(())
}
```

## Batch Workflows

```rust
use cosmolkit::{BatchErrorMode, MoleculeBatch};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let smiles = vec![
        "CCO".to_string(),
        "c1ccccc1".to_string(),
        "CC(=O)O".to_string(),
    ];

    let batch = MoleculeBatch::from_smiles_list(&smiles)
        .with_parallel_jobs(Some(8))
        .with_2d_coordinates(BatchErrorMode::Strict)?;

    let out = batch.to_smiles_list(BatchErrorMode::Strict)?;
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
the existing official InChI v1.07.5 / RDKit 2026.03.1 source port. IXA, AuxInfo
reconstruction, INCHIGEN, version-query, and extended-polymer entry points are
not exposed by this facade.

## Molecular Descriptors

The facade re-exports the source-backed descriptor functions from
`cosmolkit-core`:

```rust
use cosmolkit::{
    Molecule, calc_chi_0, calc_mol_formula, calc_mol_wt, calc_mqns,
    calc_num_aromatic_rings,
};

let molecule = Molecule::from_smiles("c1ccccc1O")?;
assert_eq!(calc_mol_formula(&molecule, false, true)?, "C6H6O");
assert!(calc_mol_wt(&molecule, false)? > 94.0);
assert_eq!(calc_num_aromatic_rings(&molecule)?, 1);
assert!(calc_chi_0(&molecule) > 0.0);
assert_eq!(calc_mqns(&molecule)?.len(), 42);
```

The documented descriptor surface includes molecular properties, connectivity
and shape indices, Lipinski and ring/stereo counts, MQN, Labute ASA, and
SlogP/SMR VSA. Supported rows and parameter combinations are checked
field-by-field against pinned RDKit golden data; unmodeled source states return
an explicit descriptor error.

### Descriptor count queries

Five read-only `Molecule` queries return RDKit-compatible count values and
are available with the `descriptors` feature:

```rust
let mol = Molecule::from_smiles("CCO")?;
assert_eq!(mol.num_heavy_atoms()?, 3);
assert_eq!(mol.total_atom_count()?, 9);
assert_eq!(mol.lipinski_hba()?, 1);
assert_eq!(mol.lipinski_hbd()?, 1);
assert_eq!(mol.fraction_csp3()?.to_bits(), 1.0_f64.to_bits());
```

`num_heavy_atoms` and `lipinski_hba` read only the topology;
`total_atom_count`, `lipinski_hbd` and `fraction_csp3` require the prepared
valence assignment cached by a sanitizing constructor and return the typed
`DescriptorReadError::MissingPreparedValence` otherwise — the queries never
create or install cache values themselves, and algorithm failures retain the
owned domain error through `Error::source`. `Molecule::num_atoms` keeps its
separate explicit-atom-row meaning; `total_atom_count` includes implicit and
explicit-property hydrogens (`includeNeighbors=false`). `lipinski_hba` is the
direct N/O count (not the general recursive `NumHBA`) and `lipinski_hbd` is
the donor-hydrogen sum on N/O (not the donor-atom count). All five run in
the parity pipeline over the 5000-record SMILES corpus under both
`remove_hs` parser policies (10000 observations per task).

## Fingerprints

The Rust facade exposes source-backed Morgan, AtomPair, Topological Torsion,
MACCS, RDKit topological, Avalon, and Layered fingerprints. ``TopologicalTorsion*`` is
the ordered atom-path torsion family; ``TopologicalFingerprint*`` remains
RDKit's distinct path/subgraph ``RDKFingerprintMol`` family. The applicable
families can also return typed provenance:

```rust
use cosmolkit::{
    AtomPairFingerprintParams, AvalonFingerprintParams, LayeredFingerprintLayers,
    LayeredFingerprintParams, Molecule,
    TopologicalFingerprintOutputRequest,
    TopologicalFingerprintParams, TopologicalTorsionFingerprintOutputRequest,
    TopologicalTorsionFingerprintParams, TopologicalTorsionFingerprintVector,
    fingerprint_topological_torsion, fingerprint_topological_torsion_with_output,
    fingerprint_topological_torsion_sparse_count,
};

let molecule = Molecule::from_smiles("c1ccccc1O")?;
let topological = molecule.fingerprint_topological(
    &TopologicalFingerprintParams::default(),
)?;
let provenance = molecule.fingerprint_topological_with_output(
    &TopologicalFingerprintParams::default(),
    TopologicalFingerprintOutputRequest {
        atom_bits: true,
        bit_info: true,
    },
)?;
let avalon = molecule.avalon_fingerprint(&AvalonFingerprintParams::default())?;
let atom_pair = molecule.fingerprint_atom_pair(&AtomPairFingerprintParams::default())?;
let layered = molecule.fingerprint_layered(&LayeredFingerprintParams {
    layers: LayeredFingerprintLayers::SUBSTRUCTURE,
    ..Default::default()
})?;
let torsion_params = TopologicalTorsionFingerprintParams::default();
let torsion_ids = fingerprint_topological_torsion_sparse_count(&molecule, &torsion_params)?;
let torsion_bits = fingerprint_topological_torsion(&molecule, &torsion_params)?;
let torsion_provenance = fingerprint_topological_torsion_with_output(
    &molecule,
    &torsion_params,
    TopologicalTorsionFingerprintOutputRequest {
        vector: TopologicalTorsionFingerprintVector::Bit,
        bit_paths: true,
        ..Default::default()
    },
)?;

assert_eq!(topological.n_bits(), 2048);
assert!(provenance.output.atom_bits.is_some());
assert_eq!(avalon.n_bits(), 512);
assert_eq!(atom_pair.n_bits(), 2048);
assert_eq!(layered.n_bits(), 2048);
assert!(!torsion_ids.nonzero_elements().is_empty());
assert_eq!(torsion_bits.n_bits(), 2048);
assert!(torsion_provenance.additional_output.is_some());
```

Topological Torsion also exposes sparse-bit and folded-count forms, ordered
``MoleculeBatch`` conveniences, shared ``FingerprintAdditionalOutput`` provenance, and
three explicitly typed legacy adapters. Invalid arguments return
``FingerprintError``; batch calculation errors retain the original record
index. Exact parity is continuously checked against pinned RDKit 2026.03.1 on
focused branch fixtures and every row of a 5,000-molecule, nine-profile
matrix. The complete ChEMBL 37 audit additionally covers all 2,897,804
mutually parseable records through 127,503,376 exact vector and provenance
comparisons. Legacy adapters preserve their historical unfolded-size and
``n_bits_per_entry`` threshold differences while delegating to the same
chemistry and vector-assembly core.

The documented topological and Avalon profiles are checked against pinned
RDKit across all 2,897,804 mutually parseable ChEMBL 37 molecules. The
full-corpus audit completed 113,014,356 exact comparisons over 14 topological
vectors, 23 Avalon vectors, and two complete topological provenance outputs
with zero mismatches. The committed 5,000-row matrices remain the continuous
regression gates for these profiles.

AtomPair is additionally checked across all 2,897,804 mutually parseable
ChEMBL 37 molecules, covering 118,809,964 comparisons over 40 vectors and one
complete provenance output per molecule with zero mismatches.

Layered exposes the six source layers, arbitrary retained source flags,
inclusive path bounds, rooted linear or branched enumeration, exact-width bit
masks, and seeded atom counts through one read-only core while preserving the
upstream ``0.7.0`` compatibility metadata. ``None`` roots mean whole-molecule
enumeration; an explicitly empty root vector enumerates no paths. Invalid
bounds, widths, count lengths, masks, and roots return ``FingerprintError``.
The complete ChEMBL 37 audit covers all 2,897,804 mutually parseable records
across 18 profiles and 52,160,472 exact comparisons with zero mismatches.
Pinned RDKit's unrooted linear branch can consume atom indices as bond indices
and terminate the process; COSMolKit deliberately uses the documented
bond-path semantics instead of reproducing that crash.

## Conformer Generation And Force Field Applications

Native conformer generation uses RDKit-aligned distance-geometry parameters.
The default value-style molecule operation uses ETKDGv3 and returns a new
molecule value. Multi-conformer generation supports deterministic seeded runs,
RMS pruning, and sequential seed expansion:

```rust
use cosmolkit::{EmbedParameters, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CC(=O)NC")?.with_hydrogens()?;

    let embedded = molecule.with_3d_conformer()?;
    println!("{}", embedded.conformers_3d().len());

    let mut params = EmbedParameters::etkdg();
    params.random_seed = 123;
    params.num_threads = 1;
    params.prune_rms_thresh = 0.5;

    let pruned = molecule.with_3d_conformers_with_params(5, params)?;
    println!("{}", pruned.conformers_3d().len());
    Ok(())
}
```

Force-field APIs operate on molecules with existing 3D conformers and return
new molecule values, so the input coordinates are left unchanged.

```rust
use cosmolkit::{
    Molecule, mmff_has_all_molecule_params, mmff_optimize_molecule,
    uff_has_all_molecule_params, uff_optimize_molecule,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let molecule = Molecule::from_smiles("CCO")?.with_hydrogens()?.sanitize()?;

    let mut builder = molecule.to_builder();
    builder.add_3d_conformer(vec![
        [0.000, 0.000, 0.000],
        [1.540, 0.000, 0.000],
        [2.100, 1.200, 0.000],
        [-0.600, 0.900, 0.000],
        [-0.600, -0.900, 0.000],
        [0.000, 0.000, 1.000],
        [1.900, -0.900, 0.000],
        [1.700, 0.000, 1.000],
        [2.900, 1.200, 0.000],
    ])?;
    let molecule = builder.build()?;

    if uff_has_all_molecule_params(&molecule)? {
        let result = uff_optimize_molecule(&molecule, 200, 10.0, -1, true)?;
        println!("UFF energy: {:.6}", result.energy);
    }

    if mmff_has_all_molecule_params(&molecule)? {
        let result = mmff_optimize_molecule(&molecule, "MMFF94", 200, 100.0, -1, true)?;
        println!("MMFF94 needs_more: {}", result.needs_more);
    }

    Ok(())
}
```

## Molecular Alignment And RMSD

Ordinary molecular alignment uses the source-backed RDKit MolAlign boundary.
Transform queries, best RMSD, coordinate-frame RMSD, and all-conformer pair
measurements are read-only. Coordinate changes are exposed only through a
value-style method or an explicit trailing-underscore method.

```rust
use cosmolkit::{AlignmentAtomMap, AlignmentParameters, Molecule};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let reference = Molecule::from_smiles("CCC")?.with_only_3d_conformer(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0]],
        true,
    )?;
    let probe = Molecule::from_smiles("CCC")?.with_only_3d_conformer(
        vec![[3.0, -2.0, 1.0], [4.0, -2.0, 1.0], [3.0, 0.0, 1.0]],
        true,
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

    let measured = probe.alignment_transform_to(&reference, &params)?;
    let (aligned, applied) = probe.with_alignment_to(&reference, &params)?;
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
