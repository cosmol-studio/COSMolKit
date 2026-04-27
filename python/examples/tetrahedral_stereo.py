"""COSMolKit tetrahedral stereo representation.

`tetrahedral_stereo_from_smiles()` exposes the internal ordered-ligand view
derived from COSMolKit's RDKit-compatible chiral tags.
"""

from cosmolkit import tetrahedral_stereo_from_smiles

stereo = tetrahedral_stereo_from_smiles("[13CH3:7][C@H](F)Cl")

for center, ligands in stereo:
    print("center:", center)
    print("ordered ligands:", ligands)
    print("reference ligand:", ligands[3])

# `None` represents an implicit hydrogen ligand.
# The first three ligands are ordered so that, when coordinates are available,
# det(ligand0 - center, ligand1 - center, ligand2 - center) is positive.
