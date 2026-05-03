"""COSMolKit usage: IO + basic properties."""

import numpy as np

from cosmolkit import Molecule

mol = Molecule.from_smiles("CCO", sanitize=True).with_2d_coords()
coords = mol.coords_2d()
print("coords shape:", coords.shape)
print("centroid:", coords.mean(axis=0))

sdf_text = mol.to_sdf_string(format="v2000")
print("SDF length:", len(sdf_text))

saved_path = mol.write_sdf_to_directory(
    "python/examples", file_name="ethanol.sdf", format="v2000"
)
print("Saved:", saved_path)

lig = Molecule.read_sdf(saved_path, sanitize=True, coordinate_dim="2d")
print("Loaded:", lig)
lig_coords = lig.coords_2d()
print("loaded coords shape:", lig_coords.shape)
print("max coordinate delta after SDF roundtrip:", np.abs(coords - lig_coords).max())
