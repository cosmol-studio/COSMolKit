"""Ideal COSMolKit usage: explicit sanitize/perception pipeline.

This script demonstrates the target public API shape only.
Current bindings are placeholders and will raise NotImplementedError.
"""

from cosmolkit import Molecule

mol = Molecule.from_smiles("c1ccccc1", sanitize=False)

report = mol.check_valence()
_ = report

mol = (
    mol.sanitize(strict=True)
    .perceive_rings()
    .perceive_aromaticity()
    .kekulize(sanitize=False)
)

_ = mol
