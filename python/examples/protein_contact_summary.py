"""Protein-chain residue and atom traversal.

Usage:
    .venv/bin/python python/examples/protein_contact_summary.py

The high-level Protein API exposes protein-only chains, residues, and atoms.
This example computes simple residue centroids and nearby residue pairs without
touching lower-level structural tables.
"""

from __future__ import annotations

from math import dist
from typing import cast

import numpy as np
import numpy.typing as npt

from cosmolkit import Protein, ProteinAtom, ProteinResidue


PDB = """\
ATOM      1  N   MET A   1      11.104  13.207   9.900  1.00 20.00           N
ATOM      2  CA  MET A   1      12.210  13.912  10.555  1.00 20.00           C
ATOM      3  C   MET A   1      13.470  13.079  10.413  1.00 20.00           C
ATOM      4  O   MET A   1      13.590  11.870  10.650  1.00 20.00           O
ATOM      5  N   GLY A   2      14.530  13.650  10.980  1.00 20.00           N
ATOM      6  CA  GLY A   2      15.790  12.920  10.910  1.00 20.00           C
ATOM      7  C   GLY A   2      16.840  13.820  11.480  1.00 20.00           C
ATOM      8  O   GLY A   2      17.980  13.450  11.520  1.00 20.00           O
ATOM      9  N   ALA B   1      10.000  15.000   8.000  1.00 20.00           N
ATOM     10  CA  ALA B   1      10.920  16.050   7.600  1.00 20.00           C
ATOM     11  C   ALA B   1      12.300  15.520   7.850  1.00 20.00           C
ATOM     12  O   ALA B   1      12.520  14.310   7.900  1.00 20.00           O
HETATM   13  O   HOH A 101      18.000  10.000   8.000  1.00 10.00           O
HETATM   14  C1  LIG C   1      18.500  11.000   8.500  1.00 10.00           C
"""


def atom_position(atom: ProteinAtom) -> tuple[float, float, float]:
    position = atom.position()
    if position is None:
        raise ValueError(f"atom {atom.name()} has no coordinates")
    return position


def residue_centroid(residue: ProteinResidue) -> npt.NDArray[np.float64]:
    coords = np.array([atom_position(atom) for atom in residue.atoms()])
    return cast(npt.NDArray[np.float64], coords.mean(axis=0))


protein = Protein.from_pdb_str(PDB)

print("models:", protein.num_models())
print("chains:", protein.num_chains())
print("residues:", protein.num_residues())
print("atoms:", protein.num_atoms())

for chain in protein.chains():
    residues = chain.residues()
    print("chain:", chain.index(), chain.kind(), "residues=", len(residues))
    for residue in residues:
        centroid = residue_centroid(residue)
        atom_names = [atom.name() for atom in residue.atoms()]
        print(
            "  residue:",
            residue.index(),
            residue.name(),
            residue.kind(),
            "atoms=",
            atom_names,
            "centroid=",
            centroid.round(3).tolist(),
        )

all_residues = protein.residues()
centroids = [residue_centroid(residue) for residue in all_residues]
print("nearby residue centroids:")
for left_index, left in enumerate(all_residues):
    for right_index in range(left_index + 1, len(all_residues)):
        distance = dist(centroids[left_index], centroids[right_index])
        if distance <= 6.0:
            right = all_residues[right_index]
            print(
                f"  {left.name()}#{left.index()} - "
                f"{right.name()}#{right.index()}: {distance:.2f} A"
            )
