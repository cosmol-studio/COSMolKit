"""COSMolKit chirality and tetrahedral stereo representation.

Atom `chiral_tag()` values keep the SMILES/RDKit-style CW and CCW path as
`ChiralTag` IntEnum values.
`Molecule.tetrahedral_stereo()` returns ordered-ligand records where the
ligand order itself carries the tetrahedral configuration. Ligand indices
match the same atom indexing used by `mol.atoms()`. Even ligand permutations
that preserve handedness are canonicalized to one numeric representative.

Specification:
https://github.com/cosmol-studio/COSMolKit/blob/main/dev/tetrahedral_stereo.md
"""

from cosmolkit import ChiralTag, Molecule

mol = Molecule.from_smiles("[13CH3:7][C@H](F)Cl")

print("isomeric smiles:", mol.to_smiles())
print("non-isomeric smiles:", mol.to_smiles(isomeric_smiles=False))

for atom in mol.atoms():
    if atom.chiral_tag() != ChiralTag.CHI_UNSPECIFIED:
        print("chiral atom:", atom.idx(), atom.chiral_tag().name)

print("chiral centers:", mol.find_chiral_centers(include_unassigned=False))

for center, ligands in mol.tetrahedral_stereo():
    print("center:", center)
    print("ordered ligands:", ligands)

opposite = Molecule.from_smiles("[13CH3:7][C@@H](F)Cl")
print("opposite ordered ligands:", opposite.tetrahedral_stereo())

with_explicit_h = mol.with_hydrogens()
print("with explicit hydrogens:", with_explicit_h.tetrahedral_stereo())

fully_substituted = Molecule.from_smiles("F[C@](Cl)(Br)I")
print("fully substituted ordered ligands:", fully_substituted.tetrahedral_stereo())
fully_substituted_opposite = Molecule.from_smiles("F[C@@](Cl)(Br)I")
print(
    "fully substituted opposite ordered ligands:",
    fully_substituted_opposite.tetrahedral_stereo(),
)

# `None` means an implicit hydrogen ligand exists but has no atom index in the
# current molecule graph. It does not mean the ligand slot is empty.
# The first three ligands are ordered so that, when coordinates are available,
# det(ligand0 - center, ligand1 - center, ligand2 - center) is positive.
