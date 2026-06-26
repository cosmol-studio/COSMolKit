"""Pickle a Molecule and restore it in memory or from disk."""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import cast

import numpy as np

from cosmolkit import Molecule


OUTPUT_DIR = Path(__file__).resolve().parent / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

mol = (
    Molecule.from_smiles("F[C@H](Cl)[13CH3:7]", sanitize=True)
    .with_hydrogens()
    .with_2d_coordinates()
)

payload = pickle.dumps(mol, protocol=pickle.HIGHEST_PROTOCOL)
restored = cast(Molecule, pickle.loads(payload))

pickle_path = OUTPUT_DIR / "molecule.pkl"
_ = pickle_path.write_bytes(payload)
restored_from_file = cast(Molecule, pickle.loads(pickle_path.read_bytes()))

print("pickle bytes:", len(payload))
print("original smiles:", mol.to_smiles(canonical=False))
print("restored smiles:", restored.to_smiles(canonical=False))
print("file restored smiles:", restored_from_file.to_smiles(canonical=False))
print("atom count:", len(restored))
print("bond count:", restored.num_bonds())
print("has 2d coordinates:", restored.has_2d_coordinates())
coordinate_delta = cast(float, np.abs(mol.coordinates_2d() - restored.coordinates_2d()).max())
print("max 2d coordinate delta:", coordinate_delta)
