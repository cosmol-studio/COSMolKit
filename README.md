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

## Installation

```bash
pip install cosmolkit
```

## Overview

**COSMolKit** is a Rust-native cheminformatics and structural biology toolkit for molecules, SMILES/SDF/MolBlock parsing, molecular graphs, conformers, coordinates, and AI-ready batch workflows.

It currently focuses on a chemistry core whose selected features are tested for RDKit-compatible behavior: SMILES parsing/writing, atom and bond feature inspection, hydrogen transforms, Kekulization, stereochemistry checks, distance-geometry bounds, SDF output, and 2D depiction.

COSMolKit is designed around ndarray-oriented structural data access, keeping molecular data efficient and natural for NumPy and PyTorch workflows.

### Copy-on-Write style transformations

COSMolKit uses deterministic, Copy-on-Write style APIs: normal molecule operations return new objects and do not mutate their inputs. This follows the modern pandas 3.0 / Polars direction of explicit dataflow instead of hidden inplace mutation. See the pandas Copy-on-Write migration guide for the related dataframe design: <https://pandas.pydata.org/docs/dev/user_guide/migration.html>.

```python
mol = Molecule.from_smiles("CCO")
mol_h = mol.with_hydrogens()

assert mol is not mol_h
```

Explicit mutation is separated into editing workflows such as `Molecule.edit()`. Internally, unchanged topology / conformer / property data can be shared so value-style transformations remain efficient.

## Current Status

COSMolKit is in early active development. The implementation is intentionally a subset and is expanded by adding well-tested behavior rather than by broad API cloning.

- **Core layout:** `cosmolkit-core` contains chemistry perception, IO, drawing, and biomolecular primitives; `cosmolkit` is the Rust facade crate; `python/` contains the PyO3 package.
- **RDKit reference:** RDKit 2026.03.1 is the active compatibility reference for selected behaviors, with `third_party/rdkit` pinned to `Release_2026_03_1` (`351f8f378f8ad6bbd517980c38896e66bf907af8`).
- **Parity coverage:** current tests cover graph features, add-H / remove-H roundtrips, tetrahedral stereo geometry, DG bounds matrices, Kekulization branches, SMILES writer branches, V2000 molblock output, and SDF V2000/V3000 roundtrips.
- **Batch-native workflows:** `MoleculeBatch` APIs support ordered molecule construction, parallel transforms, image/SDF export, and structured error handling for high-throughput datasets.
- **Python bindings:** the package exposes SMILES parsing, `Molecule.from_rdkit()`, graph/stereo inspection, value-style transforms, explicit `coords_2d()` / `coords_3d()` access, 2D/3D SDF IO for molecules with stored coordinates, DG bounds, SVG/PNG rendering, and explicit editing.
- **AI direction:** planned COSMolKit-native APIs include model-ready graph export, internal coordinates, torsion/chirality-aware diffusion helpers, and molecular tokenization. See `ai_native_features_sketch.md`.

## Python Quick Start

```python
from cosmolkit import Molecule, MoleculeBatch

mol = Molecule.from_smiles("c1ccccc1O")
mol_2d = mol.with_2d_coords()
mol_2d.write_png("phenol.png", width=400, height=300)

smiles = ["CCO", "c1ccccc1", "CC(=O)O"]
batch = MoleculeBatch.from_smiles_list(smiles, sanitize=True, errors="keep", n_jobs=8)

prepared = batch.add_hydrogens(errors="keep").compute_2d_coords(errors="keep")
prepared.to_images("molecule_images", format="png", size=(300, 300), n_jobs=8, errors="skip")
```

For more Python examples, see `python/README.md` and `python/examples/`.

## Design Principles

- **RDKit-compatible core behavior:** parity matters for molecular facts such as valence handling, aromaticity, Kekulization, stereochemistry, and file output.
- **Modern public API:** normal transformations are value-style and deterministic; explicit mutation belongs in editing workflows.
- **Rust-first performance:** heavy logic lives in Rust and is exposed to Python through PyO3 without leaking mutation ambiguity to users.
- **ndarray-oriented structural data access:** molecular data should be exposed through efficient array views that fit naturally into NumPy and PyTorch workflows.
- **Batch-native throughput:** large molecule collections should be processed as ordered batches with Rust-side parallel scheduling, minimal Python-loop overhead, and traceable per-record failures.
- **AI-ready extensions:** RDKit compatibility is the correctness floor, while COSMolKit-native graph, geometry, torsion, diffusion, and token APIs are the API ceiling.

## Roadmap

### Phase 1 — Chemistry Core Parity
**Goal:** keep the small core correct before expanding breadth

- ✅ Atom / Bond / Molecule data model
- ✅ adjacency representation
- ✅ bond order + formal charge support
- ✅ SMILES parser for the active parity corpus
- ✅ ring perception for the active parity corpus
- ✅ basic valence handling for the active parity corpus
- ✅ Kekulization for the active parity corpus
- ✅ SMILES writer parity for the active corpus
- ✅ atom and bond feature parity tests against RDKit
- ✅ explicit hydrogen expansion for the active parity corpus
- ✅ tetrahedral stereo ordered-ligand representation
- ✅ DG bounds matrix parity for the active corpus
- ✅ strict molblock V2000 coordinate/topology parity for the active corpus
- [ ] promote sanitization into a single explicit public pipeline

### Phase 2 — Chemical File I/O
**Goal:** make molecule import/export usable beyond SMILES

- [ ] MOL reader
- ✅ SDF reader with robust multi-record handling
- ✅ SDF writer with strict RDKit-compatible V2000/V3000 output
- ✅ SMILES output via RDKit-parity writer branches
- [ ] batch molecule loading
- [ ] format validation tools with precise error reporting

### Phase 2.5 — Batch-Native Processing
**Goal:** make high-throughput molecule preparation and export a core product identity

- ✅ `MoleculeBatch.from_smiles_list()` with input-order preservation
- ✅ batch transformations for sanitize, add/remove hydrogens, Kekulization, and 2D coordinates
- ✅ Rust-side parallel scheduling with configurable `n_jobs`
- ✅ structured batch errors with `errors="raise" | "keep" | "skip"`
- ✅ validity masks, error summaries, and JSON/CSV error reports
- ✅ parallel SDF and image export for large molecule collections

### Phase 3 — Python API and User Workflows
**Goal:** expose the verified Rust core through a practical Python interface

- ✅ PyO3 package scaffold
- ✅ `Molecule.from_smiles()`
- ✅ `Molecule.from_rdkit()`
- ✅ atom and bond graph access
- ✅ `Molecule.with_hydrogens()`
- ✅ `Molecule.without_hydrogens()`
- ✅ `Molecule.with_kekulized_bonds()`
- ✅ `Molecule.tetrahedral_stereo()`
- ✅ `Molecule.with_2d_coords()`
- ✅ `Molecule.to_smiles()`
- ✅ SDF read/write bindings
- ✅ SVG/PNG rendering and file export
- ✅ explicit `Molecule.edit()` workflow
- [ ] stable graph-extraction helpers for ML workflows
- [ ] explicit sanitization and error reporting API

### Phase 4 — Query, Descriptors, and Computation
**Goal:** enable practical filtering and analysis

- ✅ distance-geometry bounds matrix parity
- [ ] atom selection API
- [ ] bond selection API
- [ ] neighborhood queries
- [ ] connected component analysis
- [ ] substructure matching
- [ ] molecular formula
- [ ] molecular weight
- [ ] ring statistics
- [ ] fingerprint generation and similarity metrics

### Phase 5 — 2D Coordinates and Drawing
**Goal:** provide RDKit-drawer-like molecule depiction

- ✅ 2D coordinate generation
- ✅ SVG molecule drawer
- ✅ PNG rendering path for Python users
- ✅ embedded Noto Sans font for PNG rendering
- [ ] atom and bond annotation overlays
- [ ] stereochemistry-aware wedge/dash depiction
- [ ] visual regression tests for generated drawings

### Phase 6 — Biomolecular Structure Support
**Goal:** cover core Biopython-like structure functionality

- [ ] PDB parser
- [ ] mmCIF parser
- [ ] Structure / Model / Chain / Residue / Atom hierarchy
- [ ] alternate location handling
- [ ] insertion code handling
- [ ] HETATM parsing
- [ ] ligand extraction
- [ ] residue and chain selection utilities
- [ ] residue neighborhood queries

### Phase 7 — WASM and Browser Integration
**Goal:** enable browser-native chemistry workflows

- [ ] WASM compilation target
- [ ] JS bindings
- [ ] in-browser SMILES/SDF parsing
- [ ] lightweight molecule processing
- [ ] integration with visualization tools

### Phase 8 — AI-Native Molecular Representations
**Goal:** expose model-ready molecular data structures for modern ML workflows

- [ ] versioned `Molecule.to_graph()` export with `cosmol-v1` node and edge feature schemas
- [ ] optional graph fields for coordinates, chirality, torsions, rings, fragments, and rotatable bonds
- [ ] Python output adapters for NumPy dictionaries and graph-learning libraries
- [ ] `Molecule.to_internal_coordinates()` and `InternalCoordinates.to_cartesian()` APIs
- [ ] Z-matrix and bond-angle-torsion tree support
- [ ] ring-aware internal coordinates and torsion graph metadata
- [ ] torsion/chirality-aware diffusion utilities with periodic angle losses and sin/cos encodings
- [ ] chirality-preserving reconstruction checks and ring torsion constraints
- [ ] `Molecule.to_tokens()` with versioned graph, fragment, torsion, 3D geometry, and pharmacophore token schemes

Design sketch: `ai_native_features_sketch.md`

## Respect for RDKit

COSMolKit is developed with deep respect for RDKit and the broader open-source cheminformatics community. The goal is an independent Rust-native implementation that preserves interoperability and behavioral compatibility where appropriate, while offering a more deterministic Python API and AI-native extension surface.
