# COSMolKit Tetrahedral Stereo Representation

COSMolKit represents a tetrahedral stereocenter as `center + ordered_ligands`, instead of using only `CW/CCW`.

```rust
pub struct TetrahedralStereo {
    pub center: usize,
    pub ligands: [LigandRef; 4],
}
```

`ligands` stores the ordered four ligands around the center atom in molecule atom-index space.
If an implicit hydrogen exists, it is represented as `ImplicitH` and is usually placed at `ligands[3]`.

When 3D coordinates are available, the first three ligands should satisfy:

```text
det(pos(ligands[0]) - pos(center),
    pos(ligands[1]) - pos(center),
    pos(ligands[2]) - pos(center)) > 0
```

Thus, stereochemistry is directly defined by the ordered ligand list. `CW/CCW` is retained only as an external compatibility representation for RDKit/SMILES-style I/O.
