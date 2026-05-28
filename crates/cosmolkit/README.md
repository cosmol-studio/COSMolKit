# COSMolKit Rust

`cosmolkit` is the Rust facade crate for COSMolKit. It re-exports the molecular
model, chemistry operations, file I/O, fingerprints, drawing, batch helpers, and
protein structure APIs from `cosmolkit-core`.

## Documentation

- Rust API documentation: <https://docs.rs/cosmolkit/latest/cosmolkit/>
- Core source layout: [`../cosmolkit-core/src/README.md`](../cosmolkit-core/src/README.md)
- Python package notes: [`../../README.md`](../../README.md)
- Development policy and operation rules: [`../../dev/README.md`](../../dev/README.md)

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
its pre-call value; use the non-mutating operation when failure-preserving value
semantics are required.

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

## Examples

```bash
cargo run -p cosmolkit-core --example smiles_minimal_roundtrip
cargo run -p cosmolkit-core --example draw_svg
cargo run -p cosmolkit-core --example draw_png
cargo run -p cosmolkit-core --example sdf_to_smiles
cargo run -p cosmolkit --example protein_from_pdb
cargo run -p cosmolkit --example read_xyz
```

## Development

Core validation should use operation-contract checks:

```bash
cargo check -p cosmolkit-core --features op-contracts-strict
cargo test -p cosmolkit-core --features op-contracts-strict
cargo check -p cosmolkit-py
cargo fmt --all
```

Python binding development:

```bash
uv sync --group dev
.venv/bin/maturin develop --manifest-path python/Cargo.toml
.venv/bin/pytest
```

The facade crate should stay thin. Public Rust APIs should be exposed through
`cosmolkit` or clearly scoped public modules, while molecule mutation continues
to go through registered operations in the core.
