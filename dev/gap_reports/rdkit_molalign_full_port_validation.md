# RDKit MolAlign Full-Port Validation

## Validated Boundary

The completed boundary is ordinary RDKit `MolAlign` and
`RDNumeric::Alignments::AlignPoints` from pinned RDKit `2026.03.1`, revision
`351f8f378f8ad6bbd517980c38896e66bf907af8`. It includes:

- explicit-map and automatic first-match alignment transforms;
- best-map alignment and read-only best RMSD;
- coordinate-frame RMSD without alignment;
- all-conformer best-RMS triangular output;
- selected/all conformer alignment with ordered RMS reports;
- source conformer-ID, atom-map, weights, reflection, iteration, match-limit,
  terminal-group symmetry, and thread-count behavior;
- registered value-style and explicit trailing-underscore coordinate mutation.

O3A, MMFF/Crippen scoring, constraint O3A, and `RandomTransform` are separate
source capabilities and are not part of this boundary.

## Implementation Audit

There is one private numerical implementation at
`crates/cosmolkit-core/src/chemistry/numerics/alignment.rs`. Depiction,
distance geometry, MolAlign, and molecular coordinate transforms delegate to
its fixed-array `Transform3D`, point application, matrix multiplication,
weighted sums, covariance, quaternion, Jacobi, reflection, SSR, and translation
functions. The final audit removed the remaining caller-local 4x4 transform
helpers and added a source-tree architecture guard against their return.

The implementation is pure Rust. The pinned source call graph has no external
linear-algebra dependency, and exact covered behavior was reproduced without
an FFI backend. No `extern "C"`, backend selector, or second alignment kernel
exists in the ordinary MolAlign runtime path.

Automatic mapping reuses the sole ordinary-molecule substructure matcher with
the exact MolAlign parameter sets. Terminal conjugated-group symmetry reuses
the same clone-only helper as distance-geometry pruning. Conformer lookup
reuses one ID-based resolver. Coordinate writes are the registered
`with_alignment_to` and `with_aligned_conformers` operations under the
`coordinates.molalign` feature; read-only measurements do not enter the
mutation registry.

All copied-source markers in the MolAlign numerical, mapping, symmetry, and
conformer-selection modules are behaviorally closed for the modeled boundary.
Local complexity review found fixed-size stack matrices, source-equivalent
loop nesting and reductions, no repeated graph scan absent from the source,
and no added hot-path allocation beyond the source adapters.

## Maintained Gates

The focused oracle contains 14 source-selected rows and compares complete
transforms, RMSD, selected maps, conformer IDs/order, structured errors, and
before/after molecule state. The project-small gate covers all 152 maintained
SMILES rows. The explicit exhaustive gate covers all 5,000 maintained rows and
rotates the six public call families through the corpus. All three gates pass
against pinned RDKit with no mismatch.

Expected data is prepared reproducibly with:

```bash
.venv/bin/python tools/testdata/rdkit/generate_all.py \
  --python .venv/bin/python --suite alignment --profile smiles_5000 --jobs 4
```

The exhaustive Rust gate is:

```bash
cargo test -p cosmolkit-core --release --features op-contracts-strict \
  --test rdkit_molalign_parity \
  molalign_matches_all_5000_rows_and_rotated_operation_branches \
  -- --ignored --exact
```

## ChEMBL 37 Audit

The repository-owned `molalign` phase uses deterministic source coordinates
and rotates these call families by stable corpus row:

```text
alignment_transform
best_alignment
coordinate_rmsd
align_to
all_conformer_best_rms
align_conformers
```

The complete 128-shard execution scanned 2,897,819 ChEMBL 37 records. Fifteen
records were rejected by both parsers and 43,442 accepted records exceeded the
configured 80-atom boundary. All 2,854,362 eligible records were processed.
The phase made 11,417,207 comparisons of RMSD values, complete 4x4 transforms,
selected maps, conformer-pair order, changed coordinates, and source/reference
immutability. All 128 tasks completed with zero blocking mismatch, zero
informational mismatch, zero invalid profile, zero failed task, and no retained
finding.

The run is reproduced with:

```bash
.venv/bin/python dev/tools/chembl_parity/run.py \
  --corpus-dir /path/to/chembl37-corpus-128 \
  --run-dir /path/to/chembl37-molalign \
  --workers 96 --phase molalign
```

The corpus manifest SHA-256 is
`f8b1a516a7794c8d4428a8309c80c049b5c4f34c6c9c54381f0aebccd9ccb976` and
the corpus source SHA-256 is
`ea6181ce8dc7af41974e35b92e1febb0c9dcbe2c62f7ccc4a5d983ac19f696e7`.
Run outputs remain uncommitted; the profile, runner, comparison rules, expected
counts, source pin, and maintained regression gates are repository-owned.

## Outcome

No open behavioral gap remains inside the declared ordinary MolAlign boundary.
The explicit future boundary is O3A and its scoring/optimization dependencies;
it must not be inferred from, delegated to, or silently emulated by these APIs.
