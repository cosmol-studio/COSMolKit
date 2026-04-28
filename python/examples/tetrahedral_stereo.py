"""COSMolKit tetrahedral stereo representation.

`Molecule.tetrahedral_stereo()` exposes the internal ordered-ligand view
derived from COSMolKit's RDKit-compatible chiral tags. Ligand indices match
the same atom indexing used by `mol.atoms()`.
Specification: `tetrahedral_stereo_representation.md`.
"""

from cosmolkit import Molecule

mol = Molecule.from_smiles("[13CH3:7][C@H](F)Cl")
stereo = mol.tetrahedral_stereo()

for center, ligands in stereo:
    print("center:", center)
    print("ordered ligands:", ligands)

# `None` represents an implicit hydrogen ligand.
# The first three ligands are ordered so that, when coordinates are available,
# det(ligand0 - center, ligand1 - center, ligand2 - center) is positive.
