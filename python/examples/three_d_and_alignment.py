"""Target COSMolKit usage: 3D pipeline + pure alignment query.

This example shows planned 3D embedding, optimization, and alignment APIs that
are not implemented yet.
"""

from cosmolkit import Alignment, Molecule

mol = Molecule.from_smiles("CCO", sanitize=True)
mol3d = (
    mol.with_hydrogens()
    .embed_3d(seed=42, num_conformers=20)
    .optimize(forcefield="mmff94")
)
mol3d.write_sdf("ligand_3d.sdf")

ref = Molecule.from_smiles("CCOC(=O)N", sanitize=True)
candidates = [
    Molecule.from_smiles("CCOC(=O)N", sanitize=True),
    Molecule.from_smiles("CCN", sanitize=True),
]

best = Alignment.find_most_similar_fragment(
    reference=ref,
    candidates=candidates,
)

_ = (mol3d, best)
