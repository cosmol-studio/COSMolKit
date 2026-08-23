"""Keep complete structural data separate from protein and molecule projections."""

from cosmolkit import BioStructure


PDB = """\
ATOM      1  N   ALA A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.210  13.912  10.555  1.00 20.00           C
HETATM    3  O   HOH A   2      18.000  10.000   8.000  1.00 10.00           O
"""


structure = BioStructure.from_pdb_str(PDB)
print(structure)

# Structural writing stays on BioStructure so mixed hierarchy state is retained.
mmcif_text = structure.to_mmcif()
roundtrip = BioStructure.from_mmcif_str(mmcif_text, path="memory.cif")
print(roundtrip)

for model in structure.models():
    for chain in model.chains():
        for residue in chain.residues():
            print(residue.name(), residue.kind(), len(residue))

# Projection is explicit: the complete structure above still contains water.
protein = structure.protein()
print(protein)

# Conversion to a chemical graph is also explicit and may synthesize bonds.
molecule = structure.to_molecule(sanitize=False, remove_hs=False)
print(molecule)
