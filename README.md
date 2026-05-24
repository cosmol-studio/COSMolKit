# COSMolKit

<p align="center">
  <a href="https://github.com/cosmol-studio/COSMolKit/actions/workflows/coverage.yml">
    <img src="https://github.com/cosmol-studio/COSMolKit/actions/workflows/coverage.yml/badge.svg" alt="coverage workflow badge"/>
  </a>
  <a href="https://app.codecov.io/gh/cosmol-studio/COSMolKit">
    <img src="https://codecov.io/gh/cosmol-studio/COSMolKit/branch/main/graph/badge.svg" alt="codecov badge"/>
  </a>
  <a href="https://crates.io/crates/cosmolkit">
    <img src="https://img.shields.io/crates/v/cosmolkit.svg" alt="crates.io badge"/>
  </a>
  <a href="https://docs.rs/cosmolkit/latest/cosmolkit/">
    <img src="https://img.shields.io/docsrs/cosmolkit" alt="docs.rs badge"/>
  </a>
  <a href="https://pypi.org/project/cosmolkit/">
    <img src="https://img.shields.io/pypi/v/cosmolkit.svg" alt="pypi badge"/>
  </a>
</p>

COSMolKit is a Python molecular toolkit backed by a Rust core. It provides
value-style molecule operations, SMILES and SDF workflows, 2D depiction,
fingerprints, batch processing, and protein-focused structural biology APIs.

The library is built around explicit behavior: supported operations return
structured results, unsupported behavior fails visibly, and public molecule
transforms return new values instead of mutating their inputs.

COSMolKit is designed for array-oriented structural data access, keeping
molecular data efficient and natural for NumPy, PyTorch, and model-building
workflows.

## Documentation

- Python documentation: <https://kit.cosmol.org/>
- Rust and development notes: [`crates/cosmolkit/README.md`](crates/cosmolkit/README.md)

## Installation

```bash
pip install cosmolkit
```

## Core Concepts

- **Value-style molecules:** methods such as `with_hydrogens()`,
  `without_hydrogens()`, `with_kekulized_bonds()`, and `with_2d_coords()`
  return new molecule values.
- **Explicit errors:** invalid input and unsupported behavior are surfaced as
  errors instead of silent fallbacks.
- **Batch-native processing:** `MoleculeBatch` keeps input order, supports
  structured per-record failures, and can run batch transforms and exports with
  configurable parallelism.
- **Array-friendly data access:** coordinates, bounds matrices, fingerprints,
  and graph features are exposed in forms that fit Python numerical workflows.

### Value-Style Transformations

Normal molecule operations return new objects and do not mutate their inputs.
This follows the same explicit-dataflow direction as modern dataframe libraries:
users can reason about each transformation as a new value while the Rust core can
share unchanged internal storage efficiently.

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO")
mol_h = mol.with_hydrogens()

assert mol is not mol_h
```

## Python Quick Start

```python
from cosmolkit import Molecule, MoleculeBatch

mol = Molecule.from_smiles("c1ccccc1O")
mol_2d = mol.with_2d_coords()

print(mol_2d.to_smiles())
print(mol_2d.coords_2d())

svg = mol_2d.to_svg(width=400, height=300)
mol_2d.write_png("phenol.png", width=400, height=300)

fp = mol.fingerprint_morgan(radius=2, n_bits=2048)
print(fp.on_bits())

batch = (
    MoleculeBatch.from_smiles_list(
        ["CCO", "c1ccccc1", "CC(=O)O"],
        sanitize=True,
        errors="keep",
    )
    .with_parallel_jobs(8)
    .with_progress_bar(False)
)

prepared = batch.add_hydrogens(errors="keep").compute_2d_coords(errors="keep")
print(prepared.valid_mask())
print(prepared.to_smiles_list())

prepared.to_images(
    "molecule_images",
    format="png",
    size=(300, 300),
    errors="keep",
    filenames=["ethanol", "benzene", "acetate"],
)
```

## Protein Structures

Use `Protein` when the workflow is focused on protein chains rather than the
full structural table.

```python
from cosmolkit import Protein

protein = Protein.from_pdb("1crn.pdb")

print(protein.num_chains())
print(protein.num_residues())
print(protein.num_atoms())

for chain in protein.chains():
    print(chain.index(), chain.kind(), len(chain))
    for residue in chain.residues():
        print(residue.name(), residue.kind(), len(residue))
```

## SDF and Dataset Workflows

`SdfDataset` builds a lightweight index of SDF record byte ranges, so individual
records and chunks can be read without loading an entire file into memory.

```python
from cosmolkit import SdfDataset

dataset = SdfDataset.open("library.sdf")
print(len(dataset))

record = dataset[0]
mol = record.molecule()

for batch in dataset.batches(size=1024, errors="keep", n_jobs=8):
    smiles = batch.to_smiles_list()
```

## Feature Areas

- Molecular graph construction and inspection
- SMILES parsing and writing
- MOL/SDF reading and writing
- XYZ block reading
- Hydrogen transforms and Kekulization
- Sanitization and chemistry problem detection
- 2D coordinate generation and SVG/PNG depiction
- Morgan and Avalon fingerprints
- Distance-geometry bounds matrices
- Substructure matching and SMARTS parsing
- Ordered batch transforms and exports
- PDB/mmCIF molecule-block parsing and protein projection APIs
- Support-status metadata for public features

## Design Principles

COSMolKit aims to be Python-friendly, batch-friendly, and suitable for
model-building workflows.

- Correctness comes before breadth.
- Public transforms use value semantics.
- Mutation-capable workflows are explicit.
- Unsupported chemistry should fail clearly.
- RDKit-parity behavior is the correctness floor for supported
  cheminformatics features.
- High-throughput APIs should preserve input order and expose per-record
  failures.

## Examples

Python examples live in `python/examples/`.

## Roadmap

Status labels:

- ✅ available in the public Python API
- 🧪 implemented or partially available, still being hardened
- 🚧 planned / not yet public

### Chemistry Core

Goal: keep the supported molecular core correct before expanding breadth.

- ✅ Molecule, atom, and bond graph model
- ✅ SMILES parsing
- ✅ SMILES writing with RDKit-style writer options for supported branches
- ✅ Ring perception, valence handling, aromaticity, and Kekulization
- ✅ Hydrogen addition and removal
- ✅ Sanitization for supported chemistry workflows
- ✅ Stereochemistry inspection for supported atom and bond states
- ✅ Distance-geometry bounds matrices
- ✅ Morgan fingerprints and Tanimoto similarity
- 🧪 Avalon fingerprints
- 🧪 Substructure matching and SMARTS parsing
- 🚧 Broader descriptor APIs such as formula, molecular weight, and ring
  statistics

### File I/O and Depiction

Goal: make common molecule import, export, and visualization workflows usable
from Python.

- ✅ MOL/SDF reading
- ✅ XYZ block reading
- ✅ SDF dataset indexing for large files
- ✅ SDF writing for supported V2000/V3000 branches
- ✅ PDB block to molecule conversion
- ✅ mmCIF block to molecule conversion through the same molecule-conversion
  profile
- ✅ 2D coordinate generation
- ✅ SVG drawing
- ✅ PNG export
- 🧪 RDKit-style visual parity testing for supported depiction output
- 🚧 Annotation overlays and richer drawing customization
- 🚧 3D conformer generation and embedding APIs

### Batch-Native Workflows

Goal: make high-throughput molecule preparation and export a core product
identity.

- ✅ Ordered `MoleculeBatch.from_smiles_list()`
- ✅ Batch transforms for sanitization, hydrogens, Kekulization, and 2D
  coordinates
- ✅ Configurable parallelism with `with_parallel_jobs()`
- ✅ Configurable progress display with `with_progress_bar()`
- ✅ Per-record errors, valid masks, and error reports
- ✅ Batch SMILES, image, and SDF export paths
- 🧪 Golden parity tests for parallel batch behavior
- 🚧 More streaming and chunked dataset workflows

### Protein and Structural Biology

Goal: provide practical Biopython-like structure workflows without forcing users
through low-level structural tables.

- ✅ `Protein.from_pdb()` / `Protein.from_mmcif()` high-level entry points
- ✅ Protein chain, residue, and atom iteration
- ✅ Protein-only projection from broader structural data
- 🧪 PDB/mmCIF structural parsing
- 🚧 Selection utilities for chains, residues, atoms, and neighborhoods
- 🚧 Ligand, nucleic-acid, and mixed-structure ergonomic APIs

### Python API and ML Readiness

Goal: expose verified Rust-backed behavior through a practical Python interface.

- ✅ Value-style molecule transformations
- ✅ Graph, coordinate, fingerprint, and bounds-matrix accessors
- ✅ Python examples for drawing, SDF-to-SMILES, batch processing, and proteins
- 🧪 Type stubs and documentation coverage
- 🚧 Stable model-ready graph exports
- 🚧 NumPy / PyTorch oriented adapters
- 🚧 Molecular tokenization and AI-native geometry helpers

### Browser and Deployment

Goal: support lightweight chemistry workflows outside native Python processes.

- 🚧 WASM compilation target
- 🚧 JavaScript bindings
- 🚧 Browser-native SMILES/SDF parsing and depiction

## Respect for RDKit

COSMolKit is developed with deep respect for RDKit and the broader open-source
cheminformatics community. The goal is an independent Rust-native implementation
that preserves interoperability and RDKit-parity behavior where appropriate,
while offering a deterministic Python API and AI-native extension surface.
