"""Public Python API example: file IO and basic properties."""

import numpy as np
from cosmolkit import Molecule

KEKULE_BENZENE_MOL = """kekule_benzene
  COSMolKit      2D

  6  6  0  0  0  0  0  0  0  0999 V2000
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500   -1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  5  6  1  0
  6  1  2  0
M  END
"""

mol = Molecule.from_smiles("CCO", sanitize=True).with_2d_coordinates()
coords = mol.coordinates_2d()
print("coords shape:", coords.shape)
print("centroid:", coords.mean(axis=0))

sdf_text = mol.to_2d_sdf_string(format="v2000", include_stereo=True, kekulize=True)
print("SDF length:", len(sdf_text))

saved_path = mol.write_sdf_to_directory(
    "python/examples/output",
    file_name="ethanol.sdf",
    format="v2000",
    include_stereo=True,
    kekulize=True,
)
print("Saved:", saved_path)

lig = Molecule.read_sdf(saved_path, sanitize=True, coordinate_dim="2d")
print("Loaded:", lig)
lig_coords = lig.coordinates_2d()
print("loaded coords shape:", lig_coords.shape)
print("max coordinate delta after SDF roundtrip:", np.abs(coords - lig_coords).max())

raw_benzene = Molecule.read_mol_from_str(
    KEKULE_BENZENE_MOL,
    coordinate_dim="2d",
    sanitize=False,
)
sanitized_benzene = raw_benzene.sanitize()
print("raw MolBlock bond orders:", [bond.bond_type().name for bond in raw_benzene.bonds()])
print("delayed sanitize smiles:", sanitized_benzene.to_smiles())

explicit_h_mol = Molecule.from_smiles("CCO").with_hydrogens().with_2d_coordinates()
explicit_h_sdf = explicit_h_mol.to_2d_sdf_string(format="v2000")
kept_h = Molecule.read_sdf_from_str(
    explicit_h_sdf,
    coordinate_dim="2d",
    remove_hs=False,
)
heavy_atoms = kept_h.without_hydrogens()
print("SDF atoms before delayed hydrogen removal:", len(kept_h))
print("SDF atoms after delayed hydrogen removal:", len(heavy_atoms))
