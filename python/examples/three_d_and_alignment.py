"""Ideal COSMolKit usage: 3D pipeline + explicit alignment mutation control.

This script demonstrates the target public API shape only.
Current bindings are placeholders and will raise NotImplementedError.
"""

from cosmolkit import Alignment, Molecule

mol = Molecule.from_smiles("CCO", sanitize=True)
mol3d = (
    mol.add_hydrogens()
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
    mutate_reference=False,
    mutate_candidates=False,
)

_ = (mol3d, best)
