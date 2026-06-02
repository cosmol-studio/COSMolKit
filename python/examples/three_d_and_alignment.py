"""Public Python API example: reading 3D molfile coordinates."""

from cosmolkit import Molecule

molblock = """ethane_heavy_atoms_3d
  COSMolKit      3D

  2  1  0  0  0  0            999 V2000
   -0.7700    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.7700    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
M  END
"""

mol = Molecule.read_mol_from_str(molblock, coordinate_dim="3d")
coords = mol.coordinates_3d()

print("atoms:", mol.num_atoms())
print("bonds:", mol.num_bonds())
print("3d coords shape:", coords.shape)
print("centroid:", coords.mean(axis=0))
