"""COSMolKit usage: reading 3D SDF coordinates.

COSMolKit can preserve 3D coordinates from SDF records. 3D embedding,
optimization, and alignment are not exposed as Python APIs yet.
"""

from cosmolkit import Molecule

sdf = """methane_3d
  COSMolKit      3D

  5  4  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291    0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291   -0.6291    0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6291    0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.6291   -0.6291   -0.6291 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
  1  4  1  0
  1  5  1  0
M  END
$$$$
"""

mol = Molecule.read_sdf_from_str(sdf, coordinate_dim="3d")
coords = mol.coords_3d()

print("atoms:", mol.num_atoms())
print("3d coords shape:", coords.shape)
print("centroid:", coords.mean(axis=0))
