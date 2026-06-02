"""Public Python API example: UFF/MMFF optimization of an existing 3D conformer."""

from __future__ import annotations

import numpy as np

from cosmolkit import Molecule

ethanol_3d = """ethanol_3d
  COSMolKit      3D

  9  8  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000    0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.6000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.9000   -0.9000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    1.7000    0.0000    1.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
    2.9000    1.2000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
  1  4  1  0
  1  5  1  0
  1  6  1  0
  2  7  1  0
  2  8  1  0
  3  9  1  0
M  END
"""

mol = Molecule.read_mol_from_str(ethanol_3d, coordinate_dim="3d")
start = mol.coordinates_3d().copy()

if mol.has_uff_params():
    result = mol.with_uff_optimized(max_iters=200)
    optimized = result.molecule()

    print("UFF converged:", not result.needs_more())
    print("UFF status code:", result.status_code())
    print("UFF energy:", result.energy())
    print("coordinates changed:", not np.allclose(start, optimized.coordinates_3d()))

if mol.has_mmff_params():
    result = mol.with_mmff_optimized(max_iters=200)
    print("MMFF94 converged:", not result.needs_more())
    print("MMFF94 status code:", result.status_code())
