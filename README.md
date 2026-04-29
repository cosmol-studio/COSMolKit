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

**COSMolKit** is a Rust toolkit for core cheminformatics and structural biology data processing.

The project focuses on providing a practical and reliable subset of functionality for:

- molecular file parsing and writing  
- graph-based molecular representation  
- core chemical operations such as valence handling, aromaticity, ring perception, and Kekulization  
- structural biology file parsing such as PDB and mmCIF  

COSMolKit is designed as a **systems-level library with first-class Python bindings**, allowing it to be directly integrated into existing scientific workflows while maintaining a consistent API across Rust and Python. Leveraging Rust’s compilation model, COSMolKit can be deployed across native and WebAssembly (WASM) targets, enabling seamless use in both backend systems and browser-based environments.

A key design goal is to achieve **behavioral consistency with established tools such as RDKit** for core operations. Through extensive consistency testing and validation on real-world datasets, COSMolKit aims to reproduce identical results for critical procedures (e.g., Kekulization, aromaticity, and valence handling), making it a potential **drop-in replacement for common RDKit-based workflows**.

The overall goal is to provide a **portable, reliable, and extensible foundation** for chemistry and biomolecular structure data processing.

---

## Current Progress

COSMolKit is currently focused on the chemistry core. The repository now uses a two-crate Rust layout: `cosmolkit-core` (core implementation, including chemistry perception, IO, and biomolecular primitives) and `cosmolkit` (facade re-export crate), plus RDKit-based regression tests for SMILES parsing, atom/bond feature parity, hydrogen expansion, Kekulization, minimal MOL/SDF output, and tetrahedral stereo geometry checks.

RDKit 2026.03.1 is used as the active behavioral reference, with `third_party/rdkit` pinned to `Release_2026_03_1` (`351f8f378f8ad6bbd517980c38896e66bf907af8`). The implementation is still a subset and is being expanded by source-level parity work rather than broad API coverage.
As of 2026-04-29, the Rust workspace test suite is close to green but not fully passing against the active RDKit 2026.03.1 oracle: `cosmolkit-core` graph-feature parity, tetrahedral stereo geometry parity, kekulized topology parity, and SDF reader/writer roundtrip tests are passing, while strict V2000 molblock coordinate/topology parity still has an unresolved failure at `molblock_v2000_body_matches_rdkit_coordinates_and_topology` row 77. Tetrahedral stereo has an internal ordered-ligand representation (`TetrahedralStereo`) derived from the existing RDKit-compatible atom chiral tags, with a Rust integration test that validates positive oriented volume against RDKit ETKDGv3 coordinates (`seed=42`) on the shared chiral corpus. The representation contract is documented in `tetrahedral_stereo_representation.md`.
The repository also includes a PyO3 package under `python/`, along with a GitHub Actions workflow for building and publishing Python wheels to PyPI. The binding layer remains partial overall, but now includes graph access, `Molecule.tetrahedral_stereo()`, hydrogen add/remove APIs, and SDF read/write plus string/dir export helpers.

---

## Scope

### Chemistry core
- SMILES parsing
- SDF / MOL / MOL2 reading and writing
- atom / bond / molecule data structures
- valence and bond-order handling
- aromaticity perception
- Kekulization
- ring perception
- canonical traversal utilities
- substructure matching
- basic molecular descriptors

### Structure biology core
- PDB parsing
- mmCIF parsing
- atom / residue / chain / model hierarchy
- alternate location and insertion code handling
- hetero atom and ligand record parsing
- biomolecular structure traversal and query utilities

---

## Python Binding Strategy

Python bindings are a **first-class feature**, not an afterthought.

- All core modules are designed with binding compatibility in mind
- Exposed via PyO3
- API parity between Rust and Python is maintained where possible
- Intended to support **drop-in replacement workflows** for RDKit / Biopython in common tasks

---

## Compatibility Philosophy

A key design goal of COSMolKit is **behavioral consistency with RDKit** for core operations.

For critical functions such as:

- Kekulization  
- aromaticity perception  
- valence handling  

we aim to achieve:

> **bit-level or behavior-level equivalence with RDKit outputs**

This is ensured through:

- large-scale consistency testing against RDKit
- regression test suites on real-world molecule datasets  
- strict validation of edge cases (aromatic rings, charged systems, fused systems, etc.)

The long-term goal is to make COSMolKit usable as an **in-place replacement** for RDKit in selected workflows.

---

## Roadmap

## Phase 1 — Chemistry Core Parity
**Goal:** keep the small core correct before expanding breadth

- ✅ Atom / Bond / Molecule data model
- ✅ adjacency representation
- ✅ bond order + formal charge support
- ✅ SMILES parser for the active parity corpus
- ✅ ring perception for the active parity corpus
- ✅ basic valence handling for the active parity corpus
- ✅ Kekulization for the active parity corpus
- ✅ atom and bond feature parity tests against RDKit
- ✅ explicit hydrogen expansion for the active parity corpus
- ✅ tetrahedral stereo ordered-ligand representation
- [ ] complete the remaining strict molblock V2000 coordinate/topology parity failures
- [ ] promote sanitization into a single explicit public pipeline

Deliverable:
- parse molecules from SMILES
- run the core perception path
- match RDKit behavior for the tracked regression corpus

This phase defines the correctness baseline of the project.

---

## Phase 2 — Chemical File I/O
**Goal:** make molecule import/export usable beyond SMILES

- [ ] MOL reader
- [ ] SDF reader with robust multi-record handling
- [ ] SDF writer with strict RDKit-compatible V2000 output
- [ ] batch molecule loading
- [ ] format validation tools with precise error reporting

---

## Phase 3 — Python API and User Workflows
**Goal:** expose the verified Rust core through a practical Python interface

- ✅ PyO3 package scaffold
- ✅ `Molecule.from_smiles()`
- ✅ atom and bond graph access
- ✅ `Molecule.add_hydrogens()`
- ✅ `Molecule.remove_hydrogens()`
- ✅ `Molecule.kekulize()`
- ✅ `Molecule.tetrahedral_stereo()`
- ✅ SDF read/write bindings
- [ ] stable graph-extraction helpers for ML workflows
- [ ] explicit sanitization and error reporting API

---

## Phase 4 — Query, Descriptors, and Computation
**Goal:** enable practical filtering and analysis

- [ ] atom selection API
- [ ] bond selection API
- [ ] neighborhood queries
- [ ] connected component analysis
- [ ] substructure matching
- [ ] molecular formula
- [ ] molecular weight
- [ ] ring statistics
- [ ] fingerprint generation and similarity metrics

---

## Phase 5 — 2D Coordinates and Drawing
**Goal:** provide RDKit-drawer-like molecule depiction

- [ ] 2D coordinate generation
- [ ] SVG molecule drawer
- [ ] PNG rendering path for Python users
- [ ] atom and bond annotation overlays
- [ ] stereochemistry-aware wedge/dash depiction
- [ ] visual regression tests for generated drawings

---

## Phase 6 — Biomolecular Structure Support
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

---

## Phase 7 — WASM and Browser Integration
**Goal:** enable browser-native chemistry workflows

- [ ] WASM compilation target
- [ ] JS bindings
- [ ] in-browser SMILES/SDF parsing
- [ ] lightweight molecule processing
- [ ] integration with visualization tools

---

## Vision

COSMolKit aims to become a reliable Rust foundation for:

- chemistry file processing
- molecular graph manipulation
- core chemical perception
- biomolecular structure parsing

with:

- immediate usability via Python  
- long-term portability via WebAssembly  

In short, it is a **systems-level toolkit for chemistry and structure data**, designed for correctness, performance, and broad deployment environments.


> **Respect for RDKit**
>
> COSMolKit is developed with deep respect for RDKit and the broader open-source cheminformatics community, and aims to provide an independent Rust-native implementation while maintaining interoperability and behavioral compatibility where appropriate.
