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
  <a href="https://anaconda.org/channels/cosmol-studio/packages/cosmolkit/overview">
    <img src="https://anaconda.org/cosmol-studio/cosmolkit/badges/version.svg" alt="conda package version badge"/>
  </a>
  <a href="https://www.npmjs.com/package/@cosmol-studio/cosmolkit">
    <img src="https://img.shields.io/npm/v/@cosmol-studio/cosmolkit.svg" alt="npm badge"/>
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT license badge"/>
  </a>
</p>

[COSMolKit](https://github.com/cosmol-studio/COSMolKit) is a Rust-native cheminformatics and structural biology toolkit with first-class Python and JavaScript/TypeScript bindings. It provides molecular graph operations, SMILES/SMARTS and molecular file workflows, 2D depiction, native 3D conformer generation, UFF/MMFF optimization, fingerprints, molecular descriptors, InChI, batch processing, and protein structure APIs. Selected workflows are also available directly in the browser through the [COSMolKit Web Tools](https://tools.cosmol.org/tools).


Chemistry implementations follow pinned upstream source semantics. Reference
versions, comparison rules and completed results are recorded in
[VALIDATION.md](VALIDATION.md).

COSMolKit combines a native Rust API with Python interfaces for scientific and machine-learning workflows and JavaScript/TypeScript interfaces for browser applications. Molecular graphs, coordinates, fingerprints, bounds matrices, and structural data integrate with NumPy, PyTorch and JavaScript typed arrays.

## Features

The current chemistry reference is RDKit **2026.09.1** (Python distribution
`2026.9.1`). See the [upstream-diff audit](dev/source_provenance/rdkit_2026_03_6_to_2026_09_1.md)
for the exact revision and upstream-diff scope.

- **Modular Rust architecture:** separate crates own model values and domain
  algorithms; `cosmolkit` remains the public entry point and molecule runtime.
- **Cargo features:** select functionality such as `core`, `bio`,
  `fingerprints`, or `reaction`; prerequisites are included. Default is `full`;
  see the [Rust feature guide](crates/cosmolkit/README.md#cargo-features).
- **WebAssembly support:** run the same Rust chemistry engine in the browser
  through JavaScript bindings, without a separate algorithm implementation.
- **One API, three languages:** shared operations, defaults, results, and error
  categories, with idiomatic names for Rust, Python, and JavaScript. See
  [Public API Design](dev/public_api_design.md) for cross-language naming rules
  and shared API contracts.

Each public API is one logical entry with three projections:

| Logical name | Rust | Python | JavaScript |
|---|---|---|---|
| `from_smiles` | `Molecule::from_smiles("CCO")?` | `ck.mol_from_smiles("CCO")` | `ck.molFromSmiles("CCO")` |
| `molecular_weight` | `mol.molecular_weight()` | `mol.molecular_weight()` | `mol.molecularWeight()` |
| `to_smiles` | `mol.to_smiles()` | `mol.to_smiles()` | `mol.toSmiles()` |
| `with_hydrogens` | `mol.with_hydrogens()` | `mol.with_hydrogens()` | `mol.withHydrogens()` |
| `tpsa` | `mol.tpsa()` | `mol.tpsa()` | `mol.tpsa()` |

Learn once, then move between Python notebooks, native Rust applications, and
browser-based JavaScript without relearning the chemistry API.

## Documentation

* Python documentation: https://kit.cosmol.org/python
* JavaScript documentation: https://kit.cosmol.org/javascript
* Rust documentation: https://docs.rs/cosmolkit/latest/cosmolkit/
* Interactive Web Tools: https://tools.cosmol.org/tools

## Installation

```bash
# Rust
cargo add cosmolkit

# Python
pip install cosmolkit

# JavaScript/TypeScript (browser or bundler)
npm install @cosmol-studio/cosmolkit
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
  and graph features are exposed through NumPy and JavaScript typed arrays.
- **Source-backed 3D workflows:** conformer generation and UFF/MMFF
  optimization are available through the Python and JavaScript APIs, and atom chiral
  tags can be assigned from a selected 3D conformer with pinned-RDKit parity.

### Value-Style Transformations

Normal molecule operations return new objects and do not mutate their inputs.
This follows the same explicit-dataflow direction as modern dataframe libraries:
users can reason about each transformation as a new value while COSMolKit can
share unchanged internal storage efficiently.

```python
import cosmolkit as ck

mol = ck.mol_from_smiles("CCO")
mol_h = mol.with_hydrogens()

assert mol is not mol_h
```

## Python Quick Start

```python
import cosmolkit as ck

mol = ck.mol_from_smiles("c1ccccc1O")
mol_2d = mol.with_2d_coordinates()

print(mol_2d.to_smiles())
print(mol_2d.coordinates_2d)

mol_3d = mol.with_hydrogens().with_3d_conformer()
print(mol_3d.coordinates_3d())

svg = mol_2d.to_svg(width=400, height=300)
mol_2d.write_png("phenol.png", width=400, height=300)

fp = mol.fingerprint_morgan()
print(fp.on_bits)

array = fp.to_numpy()  # uint8, one entry per bit

batch = (
    ck.mols_from_smiles_list(
        ["CCO", "c1ccccc1", "CC(=O)O"],
    )
    .with_parallel_jobs(8)
    .with_progress_bar(False)
)

prepared = batch.with_hydrogens().with_2d_coordinates()
print(prepared.valid_mask)
print(prepared.to_smiles())

prepared.write_images("molecule_images")
```

## JavaScript Quick Start

Initialize the WASM module once. Value-style transforms leave their input
unchanged; call `free()` when WASM objects are no longer needed.

```js
import init, * as ck from "@cosmol-studio/cosmolkit";

await init();

const mol = ck.Molecule.fromSmiles("c1ccccc1O");
const mol2d = mol.with2dCoordinates();
const fp = mol.fingerprintMorgan();
const batch = ck.MoleculeBatch.fromSmilesList(["CCO", "c1ccccc1", "CC(=O)O"]);

console.log(mol.toSmiles(), mol.molecularWeight());
console.log(mol2d.toSvg(400, 300));
console.log(fp.onBits);
console.log(batch.toSmiles());

batch.free();
fp.free();
mol2d.free();
mol.free();
```

See the [JavaScript guide](https://kit.cosmol.org/javascript) for browser/CDN
initialization and more examples.

## Protein Structures

Use `BioStructure` for complete PDB/mmCIF structural data, including modeled
proteins, nucleic acids, ligands, waters, entities, models, and metadata. Use
`Protein` only when an amino-acid-only projection is intended.

```python
import cosmolkit as ck
from pathlib import Path

structure = ck.BioStructure.from_pdb(Path("complex.pdb").read_text())
print(structure.num_models, structure.num_chains, structure.num_atoms)

# Structural format conversion remains on the complete structural value.
mmcif_text = structure.to_mmcif()
roundtrip = ck.BioStructure.from_mmcif(mmcif_text)
print(roundtrip.num_atoms)
```

This example reads your `complex.pdb` file.

The protein projection is explicit and leaves the full structure available:

```python
protein = structure.protein()

print(protein.num_chains)
print(protein.num_residues)
print(protein.num_atoms)

```

## SDF and Dataset Workflows

`SdfDataset` builds a lightweight index of SDF record byte ranges, so individual
records and chunks can be read without loading an entire file into memory.
Molfile-only readers such as `Molecule.read_mol()` follow RDKit
`MolFromMolBlock` boundaries: they stop after the first `M  END` line and leave
trailing SDF data fields to the SDF APIs.

```python
import cosmolkit as ck

dataset = ck.SdfDataset.open("library.sdf")
print(len(dataset))

record = dataset[0]
mol = record.molecule

for batch in dataset.batches(size=1024, errors="keep", n_jobs=8):
    smiles = batch.to_smiles()
```

This example reads your `library.sdf` file.

## Conformer Generation And Optimization

```python
import cosmolkit as ck

mol = ck.mol_from_smiles("CC(=O)NC").with_hydrogens()

params = ck.EmbedParams(random_seed=0xF00D, num_threads=1)

embedded = mol.with_3d_conformer(params)
print(len(embedded.conformers_3d))
print(embedded.coordinates_3d())

multi = mol.with_3d_conformers(5, params)
print(len(multi.conformers_3d))

uff = embedded.uff_force_field()
print(uff.energy())
outcome = uff.minimize_(max_iterations=200)
print(uff.energy())

mmff = embedded.mmff_force_field()
print(mmff.energy())
outcome = mmff.minimize_(max_iterations=200)
print(mmff.energy())
```

`with_3d_conformer()` uses source-backed distance-geometry embedding for trusted
molecular graphs: molecules without explicit hydrogens are embedded as heavy-atom-only
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
- InChI/InChIKey conversion with structured diagnostics
- Source-backed 3D atom-chiral-tag assignment
- Typed potential-stereo analysis and lazy, source-ordered stereoisomer
  enumeration with exhaustive, bounded random, uniqueness, enhanced-group,
  and optional embedding controls
- Hydrogen transforms and Kekulization
- Sanitization and chemistry problem detection
- Source-backed Rust tautomer enumeration, canonical selection, scoring,
  callbacks, current and V1 transform catalogs, and ordered result provenance
- 2D coordinate generation and SVG/PNG depiction
- Native 3D conformer generation with DG/KDG/ETDG/ETKDG parameter presets
- Read-only molecular alignment/RMSD measurement and explicit value-style or
  trailing-underscore coordinate alignment
- UFF/MMFF optimization of generated or imported 3D conformers
- Morgan, MACCS, RDKit topological, Avalon, Pattern, AtomPair, and Topological Torsion
  fingerprints, including sparse/count
  forms, provenance, tautomer-aware Pattern hashing, 2D/3D AtomPair distances,
  and ordered batch execution
- Source-backed RDKit Layered fingerprint 0.7.0 with all six active layers,
  roots, masks, seeded atom counts, and ordered batch execution; this family
  retains upstream's experimental classification
- Distance-geometry bounds matrices
- Substructure matching and SMARTS parse metadata
- Ordered batch transforms and exports
- Python pickle round-tripping for `Molecule`
- Complete PDB/mmCIF `BioStructure` parsing, Gemmi-aligned mmCIF writing,
  protein projections, and explicit structure-to-molecule conversion
- Support-status metadata for public features

## Design Principles

COSMolKit aims to be Python- and JavaScript-friendly, batch-friendly, and
suitable for scientific and browser workflows.

- Correctness comes before breadth.
- Public transforms use value semantics.
- Mutation-capable workflows are explicit.
- **Fail-closed capability boundaries:** a separately named capability outside
  documented support returns a structured error rather than fabricated
  chemistry. This is an API design rule, not an accepted mismatch within a
  supported feature.
- RDKit-parity behavior is the correctness floor for supported
  cheminformatics features.
- High-throughput APIs should preserve input order and expose per-record
  failures.
- Reference semantics come before heuristic approximation; semantic debt is treated as a correctness risk.

## Examples

Python examples live in `python/examples/`; JavaScript examples live in
[`wasm/examples/vanilla/`](wasm/examples/vanilla/).
For the current InChI interface, see
[`python/examples/inchi_roundtrip.py`](python/examples/inchi_roundtrip.py).

## Development

See [dev/README.md](dev/README.md) for local checks and contribution rules,
[python/dev.md](python/dev.md) for Python builds, and
[wasm/dev.md](wasm/dev.md) for JavaScript/TypeScript builds.

## Acknowledgments

COSMolKit includes Rust ports of code and algorithms from:

- [RDKit](https://github.com/rdkit/rdkit)
- [Gemmi](https://github.com/project-gemmi/gemmi)

## License

COSMolKit is licensed under the [MIT License](LICENSE). Vendored sources and
externally derived test fixtures retain their upstream copyright and license
terms as documented beside those files.
