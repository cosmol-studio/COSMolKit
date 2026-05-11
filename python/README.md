# COSMolKit Python

COSMolKit is a Python package for molecule graph workflows, SMILES/SDF/MOL IO,
coordinate access, Morgan fingerprints, molecule depiction, and
high-throughput batch processing.

Current Python documentation: <https://kit.cosmol.org/>

Current note: COSMolKit already preserves supported MOL/SDF query semantics
internally in Rust, but the Python package does not yet expose a public query
AST or query-matching API. That surface is still pending design.

## API Model: Copy-On-Write (COW) Molecule Values

COSMolKit's Python API uses copy-on-write (COW) value semantics for molecule
transforms. Methods such as `with_hydrogens()`, `without_hydrogens()`,
`with_kekulized_bonds()`, and `with_2d_coords()` return a new `Molecule`; they
do not mutate the original object in place.

This is intentionally different from common RDKit Python workflows, where many
operations mutate an `RWMol` or update a molecule object directly. In COSMolKit,
keep the returned value:

```python
mol = Molecule.from_smiles("CCO")
mol_h = mol.with_hydrogens()

assert mol is not mol_h
print(mol.to_smiles())
print(mol_h.to_smiles())
```

## Installation

```bash
pip install cosmolkit
```

## Quick Start

```python
from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1O")
drawn = mol.with_2d_coords()

print(mol.to_smiles())
print(drawn.atoms()[0])

drawn.write_png("python/examples/output/phenol.png", width=400, height=300)
```

Morgan fingerprints are exposed as RDKit-style sparse bit vectors. The
``on_bits()`` output is a list of bit indexes set to 1 inside a fixed-length
binary vector, not a dense neural embedding:

```python
fp = mol.fingerprint_morgan(radius=2, n_bits=2048)

print(fp.n_bits())
print(fp.on_bits())
print(fp.tanimoto(Molecule.from_smiles("c1ccccc1").fingerprint_morgan()))
```

Additional output mirrors RDKit's Morgan provenance helpers:

```python
result = mol.fingerprint_morgan_with_output(radius=2, n_bits=2048)
info = result.additional_output()

print(result.fingerprint().on_bits())
print(info.atom_counts())
print(info.bit_info_map())
```

Chiral tags and bond orders are Python ``IntEnum`` values, so code can compare
against typed constants instead of strings:

```python
from cosmolkit import BondOrder, ChiralTag, Molecule

chiral = Molecule.from_smiles("F[C@H](Cl)Br")
print(chiral.to_smiles())
print(chiral.to_smiles(isomeric_smiles=False))

for atom in chiral.atoms():
    if atom.chiral_tag() != ChiralTag.CHI_UNSPECIFIED:
        print(atom.idx(), atom.chiral_tag().name)

for bond in Molecule.from_smiles("C=C").bonds():
    if bond.bond_type() == BondOrder.DOUBLE:
        print("double bond:", bond.begin_atom_idx(), bond.end_atom_idx())
```

Read-only maps such as `BOND_ORDER_MAP` and `CHIRAL_TAG_MAP` convert external
string names to the enum members when needed.

Use `Molecule.edit()` when you want an explicit editing workflow:

```python
editor = mol.edit()
cl = editor.add_atom("Cl")
editor.add_bond(0, cl, order="single")
mol2 = editor.commit()
```

## Batch Workflows

```python
from cosmolkit import BatchErrorMode, BatchErrorType, MoleculeBatch

smiles = ["CCO", "c1ccccc1", "not-smiles"]
batch = MoleculeBatch.from_smiles_list(
    smiles,
    errors=BatchErrorMode.KEEP,
).with_parallel_jobs(8)

for error in batch.errors():
    if error.error_type() == BatchErrorType.SMILES_PARSE:
        print(error.index(), error.message())

prepared = batch.add_hydrogens(errors=BatchErrorMode.KEEP).compute_2d_coords(
    errors=BatchErrorMode.KEEP,
)
report = prepared.to_images(
    "python/examples/output/molecule_images",
    format="png",
    errors=BatchErrorMode.SKIP,
    filenames=["ethanol", "benzene", "invalid"],
)
sdf_files = prepared.to_sdf_files(
    "python/examples/output/molecule_sdf_records",
    format="v2000",
    errors=BatchErrorMode.SKIP,
    filenames=["ethanol", "benzene", "invalid"],
)
fingerprints = prepared.fingerprint_morgan_list(n_bits=2048)

print(prepared.valid_mask())
print(prepared.errors())
print([fp.on_bits() if fp is not None else None for fp in fingerprints])
print([mol.to_smiles() if mol is not None else None for mol in prepared])
print(report)
print(sdf_files)
```

The `errors` option controls invalid records:

- `errors="raise"` raises on the first batch validation failure.
- `errors="keep"` preserves failed records and exposes structured errors.
- `errors="skip"` omits failed records from the returned result or export.

Batch SMILES output accepts formatting options:

```python
canonical = prepared.to_smiles_list(canonical=True)
explicit = prepared.to_smiles_list(
    all_bonds_explicit=True,
    all_hs_explicit=True,
)
without_maps = prepared.to_smiles_list(ignore_atom_map_numbers=True)
```

Single-molecule SMILES output accepts the same writer options:

```python
benzene = Molecule.from_smiles("c1ccccc1")
ethanol = Molecule.from_smiles("CCO")

print(benzene.to_smiles(kekule=True))
print(ethanol.to_smiles(all_bonds_explicit=True))
print(ethanol.to_smiles(canonical=False, rooted_at_atom=2))
```

`with_parallel_jobs()` configures the default worker count for later batch
operations while keeping the batch value-style. Method-level `n_jobs` can still
override it for one call, and `prepared.to_list()` returns `list[Molecule |
None]` when a Python list is more convenient than iteration.

`with_progress_bar(True)` configures Rust-side progress bars for later batch
operations. Method-level `progress_bar=True` or `False` overrides the batch
default for a single call, and the progress output is written to stderr.

Directory exports accept optional `filenames` lists. Entries are aligned with
the batch records, `None` keeps the default numbered filename, and missing
extensions are filled from the selected output format.

## SDF and Arrays

```python
mol = Molecule.from_smiles("CCO").with_2d_coords()

sdf_text = mol.to_2d_sdf_string(format="v2000", include_stereo=True, kekulize=True)
restored = Molecule.read_sdf_from_str(sdf_text, coordinate_dim="2d")

coords = restored.coords_2d()
bounds = restored.dg_bounds_matrix()

print(coords.shape)
print(bounds.shape)
```

For multi-record SDF files, choose the reader by workload:

```python
from cosmolkit import MoleculeBatch, SdfDataset, SdfReader

# Convenience: load the whole file into one in-memory batch.
batch = MoleculeBatch.read_sdf("input.sdf", errors="keep", progress_bar=True)

# Large seekable files: build a lightweight index, then read records on demand.
dataset = SdfDataset.open("large.sdf")
record = dataset[12345]
head = dataset[:1024]
for batch in dataset.batches(size=1024, errors="keep", progress_bar=True):
    ...

# One-pass stream-style processing.
for batch in SdfReader.open("large.sdf").batches(size=1024, errors="keep"):
    ...
```

Use `MoleculeBatch.read_sdf()` when you want all molecules in memory. Use
`SdfDataset` when you need random access, accurate record-count progress, or
bounded-memory batches. `dataset[i]` returns one `SdfRecord`; slices and index
lists return `MoleculeBatch`. Use `SdfReader` for simple forward-only pipelines.

## Main Features

- SMILES parsing and writing with `Molecule.from_smiles()` and `to_smiles()`
- RDKit-style SMILES writer options on single molecules and batches
- copy-on-write molecule value semantics for transforms
- SDF and MOL file/string IO
- atom and bond feature inspection
- hydrogen add/remove transforms
- Kekule bond representation
- CW/CCW chiral tags, chiral centers, and tetrahedral stereo inspection
- 2D coordinate generation and NumPy coordinate arrays
- distance-geometry bounds matrix export
- Morgan fingerprint bit vectors, Tanimoto similarity, and AdditionalOutput
- SVG and PNG molecule depictions
- ordered batch construction, transformation, filtering, and export
- batch Morgan fingerprint generation with Rust-side parallel scheduling
- batch-level default parallelism with `with_parallel_jobs()`
- Python-style batch iteration, slicing, integer-list selection, and boolean-mask selection
- custom per-record filenames for batch image and SDF directory exports
- explicit molecule editing with `Molecule.edit()`
