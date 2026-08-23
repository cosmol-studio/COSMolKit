# COSMolKit

`cosmolkit` is the public Rust API for COSMolKit, a Rust-native cheminformatics and structural biology toolkit. It provides molecular graphs, SMILES/SMARTS processing, molecular file IO, fingerprints, descriptors, 2D depiction, native 3D conformer generation, UFF/MMFF optimization, InChI, substructure search, batch workflows, and protein structure APIs. The crate is a lightweight facade over `cosmolkit-core` and related COSMolKit components, providing the primary Rust import surface without hiding the underlying modules.

## Documentation

- Rust API documentation: <https://docs.rs/cosmolkit/latest/cosmolkit/>
- Core source layout: [`../cosmolkit-core/src/README.md`](../cosmolkit-core/src/README.md)
- Python package notes: [`../../README.md`](../../README.md)

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

The Rust facade exposes the four audited scalar InChI APIs directly:

```rust
use cosmolkit::{Molecule, inchi_to_inchi_key, mol_from_inchi, mol_to_inchi};

let molecule = Molecule::from_smiles("C")?;
let generated = mol_to_inchi(&molecule, None)?;
assert_eq!(generated.inchi, b"InChI=1S/CH4/h1H4");

let key = inchi_to_inchi_key(&generated.inchi)?;
assert_eq!(key.key, b"VNWKTOKETHGBQD-UHFFFAOYSA-N");

let parsed = mol_from_inchi(&generated.inchi, false, false)?;
assert!(parsed.molecule.is_some());
```

Pinned official InChI v1.07.5 and RDKit 2026.03.1 establish exact parity for
source-defined behavior in this boundary. Official-C undefined behavior on the
audited `NormalizeAndCompare` initial-allocation path is mapped to a
deterministic structured allocation error. MolBlock, SDF/V3000, IXA, AuxInfo,
INCHIGEN, version-query, and extended-polymer InChI APIs are not exposed.

## Molecular Descriptors

The facade re-exports the source-backed descriptor functions from
`cosmolkit-core`:

```rust
use cosmolkit::{Molecule, calc_mol_formula, calc_mol_wt, calc_num_aromatic_rings};

let molecule = Molecule::from_smiles("c1ccccc1O")?;
assert_eq!(calc_mol_formula(&molecule, false, true)?, "C6H6O");
assert!(calc_mol_wt(&molecule, false)? > 94.0);
assert_eq!(calc_num_aromatic_rings(&molecule)?, 1);
```

The documented descriptor surface is stable. Supported rows and parameter
combinations are checked field-by-field against pinned RDKit golden data;
unmodeled source states return an explicit descriptor error.

## Fingerprints

The Rust facade exposes source-backed Morgan, MACCS, RDKit topological, Avalon,
and AtomPair fingerprints. Topological and AtomPair fingerprints can also
return typed provenance:

```rust
use cosmolkit::{
    AtomPairFingerprintParams, AvalonFingerprintParams, Molecule,
    TopologicalFingerprintOutputRequest, TopologicalFingerprintParams,
};

let molecule = Molecule::from_smiles("c1ccccc1O")?;
let topological = molecule.topological_fingerprint(
    &TopologicalFingerprintParams::default(),
)?;
let provenance = molecule.topological_fingerprint_with_output(
    &TopologicalFingerprintParams::default(),
    TopologicalFingerprintOutputRequest {
        atom_bits: true,
        bit_info: true,
    },
)?;
let avalon = molecule.avalon_fingerprint(&AvalonFingerprintParams::default())?;
let atom_pair = molecule.atom_pair_fingerprint(&AtomPairFingerprintParams::default())?;

assert_eq!(topological.n_bits(), 2048);
assert!(provenance.output.atom_bits.is_some());
assert_eq!(avalon.n_bits(), 512);
assert_eq!(atom_pair.n_bits(), 2048);
```

The documented topological, Avalon, and AtomPair profiles are checked against
pinned RDKit across all 2,897,804 mutually parseable ChEMBL 37 molecules. The
topological/Avalon audit completed 113,014,356 exact comparisons; the AtomPair
audit completed another 118,809,964 comparisons over 40 vectors and one
complete provenance output per molecule. Both completed with zero mismatches.
The committed 5,000-row matrices remain the continuous regression gates for
these profiles.

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

## Examples

```bash
cargo run -p cosmolkit --example smiles_write_options
cargo run -p cosmolkit --example draw_svg
cargo run -p cosmolkit --example draw_png
cargo run -p cosmolkit --example sdf_to_smiles
cargo run -p cosmolkit --example protein_from_pdb
cargo run -p cosmolkit --example read_xyz
cargo run -p cosmolkit --example conformer_generation
cargo run -p cosmolkit --example forcefield_optimization
```

## Development

Core validation should use operation-contract checks:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --release --features op-contracts-strict
cargo check -p cosmolkit-py
cargo fmt --all
```

Use debug-profile test filters for small local iterations. Use release mode with
`op-contracts-strict` for large local runs, parity suites, and CI; release-mode
testing keeps operation contracts and runtime invariants enabled through the
strict feature set.

Python binding development:

```bash
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
.venv/bin/pytest
```

The facade crate should stay thin. Public Rust APIs should be exposed through
`cosmolkit` or clearly scoped public modules, while molecule mutation continues
to go through registered operations in the core.
