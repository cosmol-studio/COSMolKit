"""Target COSMolKit sanitize/perception pipeline.

This example shows the intended explicit workflow shape. The sanitize and
perception APIs used here are not implemented yet.
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
