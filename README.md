# COSMolKit

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

RDKit 2025.03.5 is used as the active behavioral reference. The implementation is still a subset and is being expanded by source-level parity work rather than broad API coverage.
As of 2026-04-28, `cosmolkit-core` RDKit graph-feature parity tests are passing for direct and explicit-hydrogen molecules. Tetrahedral stereo now also has an internal ordered-ligand representation (`TetrahedralStereo`) derived from the existing RDKit-compatible atom chiral tags, with a Rust integration test that validates positive oriented volume against RDKit ETKDGv3 coordinates (`seed=42`) on the shared chiral corpus. The representation contract is documented in `tetrahedral_stereo_representation.md`. Molblock parity remains in progress: kekulized topology parity is currently passing on the shared corpus, while strict V2000 coordinate parity currently fails first at row 52 (`COC(=O)c1ccc(NC(=O)c2cccc(C)c2)cc1`).
The repository also includes a PyO3 package under `python/`, along with a GitHub Actions workflow for building and publishing Python wheels to PyPI. The binding layer is still partial, but now exposes lower-level molecule graph access and `Molecule.tetrahedral_stereo()` for the new stereo representation.

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

## Phase 0 — Project Foundation
**Goal:** establish the minimal architecture and internal data model

- ✅ Rust workspace layout
- [ ] core error types
- [ ] shared chemical element table
- [ ] common indexing and identifier types
- [ ] serialization strategy
- ✅ test corpus for molecule parity work
- ✅ PyO3 binding scaffold

---

## Phase 1 — Minimal Usable Core (Critical Start)
**Goal:** achieve a small but complete and verifiable workflow

- [ ] SDF reader (robust, streaming-capable)
- ✅ Atom / Bond / Molecule struct
- ✅ adjacency representation
- ✅ bond order + formal charge support
- [ ] basic valence handling (in progress for tested subset)
- [ ] Kekulization (in progress for tested subset)

Deliverable:
- load molecules from SDF
- perform Kekulization
- return identical results (validated) to RDKit

This phase defines the **correctness baseline** of the project.

---

## Phase 2 — Chemical Perception Core
**Goal:** complete the minimal RDKit-equivalent perception pipeline

- [ ] implicit hydrogen assignment (in progress for tested subset)
- [ ] aromaticity detection
- [ ] ring perception (in progress for tested subset)
- [ ] sanitization pipeline
- [ ] consistency testing vs RDKit (in progress: `cosmolkit-core` graph-feature parity and kekulized molblock topology parity are green; tetrahedral stereo ETKDG geometry check added; strict V2000 coordinate parity still being closed)

Deliverable:
- molecules can be fully sanitized and normalized
- results match RDKit for standard datasets

---

## Phase 3 — Chemical File I/O Expansion
**Goal:** broaden usable chemistry workflows

- [ ] SMILES parser (in progress for tested subset)
- [ ] MOL reader
- [ ] SDF writer (minimal V2000 output in progress; strict coordinate parity still open)
- [ ] tetrahedral stereo API promotion from derived view to first-class internal geometry/stereo primitive
- [ ] batch molecule loading
- [ ] format validation tools

---

## Phase 4 — Query and Basic Computation
**Goal:** enable practical filtering and analysis

- [ ] atom selection API
- [ ] bond selection API
- [ ] neighborhood queries
- [ ] connected component analysis
- [ ] substructure matching (initial version)
- [ ] molecular formula
- [ ] molecular weight
- [ ] ring statistics

---

## Phase 5 — Biomolecular Structure Parsing
**Goal:** cover core Biopython-like structure functionality

- [ ] PDB parser
- [ ] mmCIF parser
- [ ] Structure / Model / Chain / Residue / Atom hierarchy
- [ ] alternate location handling
- [ ] insertion code handling
- [ ] HETATM parsing
- [ ] coordinate extraction utilities

---

## Phase 6 — Biomolecular Utilities
**Goal:** support real preprocessing workflows

- [ ] chain filtering
- [ ] residue selection
- [ ] ligand extraction
- [ ] water / ion filtering
- [ ] coordinate range utilities
- [ ] residue neighborhood queries

---

## Phase 7 — Stability and Ecosystem
**Goal:** make COSMolKit production-ready

- [ ] Python API stabilization
- [ ] documentation and examples
- [ ] benchmark suite
- [ ] dataset-scale performance testing
- [ ] CLI tools

---

## Phase 8 — WASM Support (Forward-Looking)
**Goal:** enable browser-native scientific computation

- [ ] WASM compilation target
- [ ] JS bindings
- [ ] in-browser SDF parsing
- [ ] lightweight molecule processing
- [ ] integration with visualization tools

---

## Non-Goals (Early Development)

- conformer generation
- docking
- scoring
- virtual screening
- molecular simulation
- protein–ligand interaction modeling
- deep learning model integration

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
