"""COSMolKit chirality and tetrahedral stereo representation.

Atom `chiral_tag()` values keep the SMILES/RDKit-style CW and CCW path.
`Molecule.tetrahedral_stereo()` is an additional ordered-ligand view derived
from those chiral tags. Ligand indices match the same atom indexing used by
`mol.atoms()`.
Specification: `tetrahedral_stereo_representation.md`.
"""

from cosmolkit import Molecule

mol = Molecule.from_smiles("[13CH3:7][C@H](F)Cl")

print("isomeric smiles:", mol.to_smiles())
print("non-isomeric smiles:", mol.to_smiles(isomeric_smiles=False))

for atom in mol.atoms():
    if atom.chiral_tag() != "CHI_UNSPECIFIED":
        print("chiral atom:", atom.idx(), atom.chiral_tag())

print("chiral centers:", mol.find_chiral_centers(include_unassigned=False))

for center, ligands in mol.tetrahedral_stereo():
    print("center:", center)
    print("ordered ligands:", ligands)

# `None` represents an implicit hydrogen ligand.
# The first three ligands are ordered so that, when coordinates are available,
# det(ligand0 - center, ligand1 - center, ligand2 - center) is positive.
