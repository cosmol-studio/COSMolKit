"""COSMolKit usage: ndarray-oriented structural data access."""

import numpy as np

from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1O").with_2d_coords()

coords = mol.coords_2d()
bounds = mol.dg_bounds_matrix()

print("coords type:", type(coords).__name__)
print("coords shape:", coords.shape)
print("centroid:", coords.mean(axis=0))
print("centered coords:", coords - coords.mean(axis=0))

print("DG bounds type:", type(bounds).__name__)
print("DG bounds shape:", bounds.shape)
print("d(0,1):", bounds[0, 1])
print("upper bounds row sums:", bounds.sum(axis=1))
print("diagonal is zero:", np.allclose(np.diag(bounds), 0.0))
