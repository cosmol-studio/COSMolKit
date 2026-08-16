# COSMolKit — Rust-native cheminformatics toolkit

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
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT license badge"/>
  </a>
</p>

[COSMolKit](https://github.com/cosmol-studio/COSMolKit) is a Rust-native cheminformatics and structural biology toolkit with first-class Python bindings. It provides molecular graph operations, SMILES/SMARTS and molecular file workflows, 2D depiction, native 3D conformer generation, UFF/MMFF optimization, fingerprints, molecular descriptors, InChI, batch processing, and protein structure APIs.

For supported cheminformatics operations, RDKit-compatible behavior is treated as the correctness floor. COSMolKit uses boundary-scoped parity claims: a feature is considered parity-covered only when its documented reference surface passes the required exact or numerical comparisons. Fixed reference oracles, source-backed implementations where reference semantics require them, committed regression corpora, and explicit unsupported-state handling are used together; aggregate success rates or approximate similarity are not treated as substitutes for behavioral parity.

COSMolKit combines a native Rust API with Python interfaces designed for array-oriented scientific and machine-learning workflows. Molecular graphs, coordinates, fingerprints, bounds matrices, and structural data are exposed in forms suitable for NumPy, PyTorch, dataset processing, and model-building pipelines.

## Documentation

- Python documentation: <https://kit.cosmol.org/>
- Rust crate notes: [`crates/cosmolkit/README.md`](crates/cosmolkit/README.md)
- Feature parity scope: [`dev/parity_scope.md`](dev/parity_scope.md)

## Validation Status

COSMolKit is extensively agent-driven, but compatibility-critical chemistry is not independently reimplemented from descriptions or approximated heuristically. Where reference behavior matters, implementations follow source-backed ports with explicit operation contracts, pinned RDKit/upstream oracles, and strict parity gates. Agents accelerate implementation; the reference source and parity tests constrain what counts as correct. COSMolKit does not infer reference semantics from final-output agreement: it source-ports reference behavior, verifies intermediate molecular states and invariants, and makes mutation explicit—because matching outputs alone cannot make a heuristic reimplementation trustworthy.

A parity claim is made only for a documented behavioral boundary that closes exactly or within an explicit numerical tolerance. **99% or 99.9% agreement remains unfinished when mismatches exist.** Unsupported branches fail explicitly rather than returning plausible substitutes. This source-first workflow is designed to minimize **semantic debt**: small heuristic deviations are not allowed to accumulate into undocumented behavioral differences that later require layers of corrective patches. By resolving reference semantics at implementation time and locking them with reproducible parity tests, COSMolKit aims to prevent locally reasonable approximations from becoming systemic chemistry errors discovered only through downstream corpus testing.

Current parity-covered surfaces include molecular graph state, SMILES branch behavior, distance geometry, ETKDG, UFF/MMFF, molecular descriptors, Morgan/MACCS fingerprints, molecular I/O, depiction, and the supported InChI APIs. Representative gates include:

* **ETKDG:** bounds matrices to `1e-8`; covered fixed-seed coordinates to `1e-6` versus pinned RDKit.
* **UFF/MMFF:** 5,000 molecules for parameter availability, MMFF atom types and charges; curated cases additionally require energy, every gradient component, optimizer result, and final coordinates within `1e-6`.
* **Descriptors:** covered floating-point outputs are checked against pinned RDKit bit patterns.
* **Fingerprints:** supported Morgan and MACCS branches require exact RDKit equality; unfinished APIs do not substitute approximate vectors.
* **InChI:** the production implementation is pure Rust and source-backed against official InChI v1.07.5 semantics; the strict 5,000-molecule profile requires exact InChI equality on every row. Defined upstream behavior is reproduced; upstream undefined behavior is isolated and mapped to explicit structured Rust errors rather than disguised as parity.

See [`dev/parity_scope.md`](dev/parity_scope.md) for the exact feature boundaries, reference versions, comparison fields, tolerances, and unfinished surfaces.

## Installation

```bash
pip install cosmolkit
```

## Core Concepts

- **Value-style molecules:** methods such as `with_hydrogens()`,
  `without_hydrogens()`, `with_kekulized_bonds()`, and `with_2d_coordinates()` return new molecule values, keeping topology-changing operations explicit and preventing derived chemistry state from being silently invalidated.
- **Explicit mutation:** in-place `Molecule` operations always end with `_`.
  The trailing underscore has no other public `Molecule` meaning.
- **Explicit errors:** invalid input and unsupported behavior are surfaced as
  errors instead of silent fallbacks.
- **Batch-native processing:** `MoleculeBatch` keeps input order, supports
  structured per-record failures, and can run batch transforms and exports with
  configurable parallelism.
- **Array-friendly data access:** coordinates, bounds matrices, fingerprints,
  and graph features are exposed in forms that fit Python numerical workflows.
- **Source-backed 3D workflows:** conformer generation and UFF/MMFF
  optimization are available through the public Python API, and atom chiral
  tags can be assigned from a selected 3D conformer with pinned-RDKit parity.

### Value-Style Transformations

Normal molecule operations return new objects and do not mutate their inputs.
This follows the same explicit-dataflow direction as modern dataframe libraries:
users can reason about each transformation as a new value while COSMolKit can
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
mol_2d = mol.with_2d_coordinates()

print(mol_2d.to_smiles())
print(mol_2d.coordinates_2d())

mol_3d = mol.with_hydrogens().with_3d_conformer()
print(mol_3d.coordinates_3d().shape)

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

prepared = batch.with_hydrogens(errors="keep").with_2d_coordinates(errors="keep")
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
Molfile-only readers such as `Molecule.read_mol()` follow RDKit
`MolFromMolBlock` boundaries: they stop after the first `M  END` line and leave
trailing SDF data fields to the SDF APIs.

```python
from cosmolkit import SdfDataset

dataset = SdfDataset.open("library.sdf")
print(len(dataset))

record = dataset[0]
mol = record.molecule()

for batch in dataset.batches(size=1024, errors="keep", n_jobs=8):
    smiles = batch.to_smiles_list()
```

## Conformer Generation And Optimization

```python
from cosmolkit import EmbedParameters, Molecule

mol = Molecule.from_smiles("CC(=O)NC").with_hydrogens()

params = EmbedParameters.etkdg_v3()
params.random_seed = 0xF00D
params.num_threads = 1
params.track_failures = True

embedded = mol.with_3d_conformer(params)
print(embedded.num_conformers())
print(embedded.coordinates_3d().shape)
print(params.failures)

multi = mol.with_3d_conformers(5, params)
print(multi.num_conformers())

if embedded.has_uff_params():
    uff = embedded.with_uff_optimized(max_iters=200)
    print(uff.energy())

if embedded.has_mmff_params():
    mmff = embedded.with_mmff_optimized(max_iters=200)
    print(mmff.needs_more())
```

`with_3d_conformer()` follows RDKit's ETKDG behavior for trusted molecular
graphs: molecules without explicit hydrogens are embedded as heavy-atom-only
conformers instead of failing or automatically adding hydrogens. Calling
`with_hydrogens()` first is recommended for all-atom geometry, force-field
optimization, and hydrogen-bond-sensitive workflows. Coordinate-only inputs
such as XYZ blocks do not contain a bond topology and are not valid ETKDG
inputs until a trusted graph has been constructed.

## Feature Areas

- Molecular graph construction and inspection
- SMILES parsing and writing
- MOL/SDF reading and writing
- MOL2 reading with RDKit-style `Mol2ParserParams`
- XYZ block reading
- Four scalar InChI APIs with exact source-defined official-C/RDKit parity and structured errors
- Stable 3D atom-chiral-tag assignment with exact pinned-RDKit full-state parity
- Hydrogen transforms and Kekulization
- Sanitization and chemistry problem detection
- 2D coordinate generation and SVG/PNG depiction
- Native 3D conformer generation with DG/KDG/ETDG/ETKDG parameter presets
- UFF/MMFF optimization of generated or imported 3D conformers
- Morgan and MACCS fingerprints for the validated exact-parity branches
- Distance-geometry bounds matrices
- Substructure matching and SMARTS parse metadata
- Ordered batch transforms and exports
- Python pickle round-tripping for `Molecule`
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
- Reference semantics come before heuristic approximation; semantic debt is treated as a correctness risk.

## Examples

Python examples live in `python/examples/`.
For the current InChI interface, see
[`python/examples/inchi_roundtrip.py`](python/examples/inchi_roundtrip.py).

## Development

Small focused Rust test filters may use the default debug profile while
iterating:

```bash
cargo test -p cosmolkit-core --features op-contracts-strict <test-filter>
```

Large local runs, parity suites, and CI tests should use release mode with the
same strict feature set:

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict
```

Release-mode testing keeps operation contracts and runtime invariants enabled
through `op-contracts-strict`; optimized release builds for distribution use
default features unless explicit runtime checks are requested.

## Roadmap

Status labels:

- ✅ stable public functionality within its documented supported scope
- 🧪 public experimental feature; available, but its behavior or API may change
- 🚧 planned or not yet public

The ✅ status applies to the documented COSMolKit scope. It does not claim that
every API or input branch from an upstream reference library is implemented;
behavior outside that scope must continue to fail explicitly.

### Chemistry Core

Goal: keep the supported molecular core correct before expanding breadth.

- ✅ Molecule, atom, and bond graph model
- ✅ SMILES parsing
- ✅ SMILES writing with RDKit-style writer options for supported branches
- ✅ Ring perception, valence handling, aromaticity, and Kekulization
- ✅ Hydrogen addition and removal
- ✅ Sanitization for supported chemistry workflows
- ✅ Stereochemistry inspection for supported atom and bond states
- ✅ Atom chiral-tag assignment from selected 3D conformers, with exact
  pinned-RDKit `assignChiralTypesFrom3D` parity across 77 fixed full-state
  oracle records
- ✅ Distance-geometry bounds matrices
- ✅ Native 3D conformer generation and UFF/MMFF post-optimization for
  supported molecules
- ✅ `Molecule.to_inchi()`, `Molecule.to_inchi_key()`, `inchi_to_key()`, and
  `Molecule.from_inchi()` for source-defined behavior; official-C undefined
  allocation behavior returns a structured error
- ✅ Morgan fingerprints and Tanimoto similarity for the validated exact-parity branches
- ✅ MACCS fingerprints for the validated exact raw/public projection
- 🚧 RDKFingerprint/topological and Avalon fingerprints fail closed until
  the [source-exact follow-up plan](dev/plans/rdkit_topological_avalon_fingerprint_port_plan.md)
  is completed with exact-bit parity
- ✅ Substructure matching and Python SMARTS parse metadata
- ✅ Molecular descriptors: average/exact molecular weight, formula, H-bond
  donor/acceptor counts, fraction Csp3, Crippen logP/MR, TPSA, aromatic-ring
  count, rotatable-bond modes, and QED for the documented parameter space

### File I/O and Depiction

Goal: make common molecule import, export, and visualization workflows usable
from Python.

- ✅ MOL/SDF reading
- ✅ MOL2 reading
- ✅ XYZ block reading
- ✅ SDF dataset indexing for large files
- ✅ SDF writing for supported V2000/V3000 branches
- ✅ PDB block to molecule conversion
- ✅ mmCIF block to molecule conversion through the same molecule-conversion
  profile
- ✅ 2D coordinate generation
- ✅ SVG drawing
- ✅ PNG export
- ✅ RDKit-style visual parity testing for supported depiction output
- 🚧 Annotation overlays and richer drawing customization
- ✅ 3D conformer generation and embedding APIs

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
- ✅ Golden parity tests for parallel batch behavior
- 🚧 More streaming and chunked dataset workflows

### Protein and Structural Biology

Goal: provide practical Biopython-like structure workflows without forcing users
through low-level structural tables.

- ✅ `Protein.from_pdb()` / `Protein.from_mmcif()` high-level entry points
- ✅ Protein chain, residue, and atom iteration
- ✅ Protein-only projection from broader structural data
- ✅ PDB/mmCIF structural parsing
- 🚧 Selection utilities for chains, residues, atoms, and neighborhoods
- 🚧 Ligand, nucleic-acid, and mixed-structure ergonomic APIs

### Python API and ML Readiness

Goal: expose verified molecular behavior through a practical Python interface.

- ✅ Stable value-style mutation contract for public molecule transformations
- ✅ Graph, coordinate, fingerprint, descriptor, and bounds-matrix accessors
- ✅ Python examples for drawing, SDF-to-SMILES, pickle round-tripping, batch
  processing, and proteins
- ✅ Type stubs and documentation coverage
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
cheminformatics community. The goal is a Rust-native implementation that preserves interoperability and faithfully ports reference behavior where appropriate, while offering a deterministic Python API and AI-native extension surface.

## License

COSMolKit is licensed under the [MIT License](LICENSE). Vendored sources and
externally derived test fixtures retain their upstream copyright and license
terms as documented beside those files.
