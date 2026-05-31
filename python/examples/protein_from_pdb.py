"""Public Python API example: protein subset traversal.

Usage:
    .venv/bin/python python/examples/protein_from_pdb.py
"""

from cosmolkit import Protein


PDB = """\
ATOM      1  N   MET A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  MET A   1      12.210  13.912  10.555  1.00 20.00           C
ATOM      3  C   MET A   1      13.470  13.079  10.413  1.00 20.00           C
ATOM      4  N   GLY A   2      14.530  13.650  10.980  1.00 20.00           N
ATOM      5  CA  GLY A   2      15.790  12.920  10.910  1.00 20.00           C
HETATM    6  O   HOH A   3      18.000  10.000   8.000  1.00 10.00           O
HETATM    7  C1  LIG B   1      18.500  11.000   8.500  1.00 10.00           C
"""


protein = Protein.from_pdb_str(PDB)
print("chains:", protein.num_chains())
print("residues:", protein.num_residues())
print("atoms:", protein.num_atoms())

first_chain = protein[0]
residues = [residue for residue in first_chain.residues()]
print("residue names:", [residue.name() for residue in residues])

for residue in residues:
    print(residue.name(), [atom.name() for atom in residue.atoms()])
