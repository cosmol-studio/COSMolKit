#!/usr/bin/env python3
"""Generate RDKit XYZ reader golden data from SMILES-derived XYZ blocks."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger

EXPECTED_RDKIT_VERSION = "2026.3.1"


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield line


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "xyz_block": None,
            "atomic_numbers": None,
            "coordinates": None,
            "comment": None,
            "num_bonds": None,
            "error": "MolFromSmiles returned None",
        }

    mol = Chem.AddHs(mol, addCoords=False)
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        atomic_num = atom.GetAtomicNum()
        conformer.SetAtomPosition(
            idx,
            (
                float(idx),
                float((idx % 7) - 3) * 0.25,
                float((atomic_num % 11) - 5) * 0.125,
            ),
        )
    mol.AddConformer(conformer, assignId=True)

    xyz_block = Chem.MolToXYZBlock(mol)
    xyz_mol = Chem.MolFromXYZBlock(xyz_block)
    if xyz_mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "xyz_block": xyz_block,
            "atomic_numbers": None,
            "coordinates": None,
            "comment": None,
            "num_bonds": None,
            "error": "MolFromXYZBlock returned None",
        }

    conf = xyz_mol.GetConformer()
    comment = (
        xyz_mol.GetProp("_FileComments")
        if xyz_mol.HasProp("_FileComments")
        else None
    )
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "xyz_block": xyz_block,
        "atomic_numbers": [atom.GetAtomicNum() for atom in xyz_mol.GetAtoms()],
        "coordinates": [
            [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
            for i in range(xyz_mol.GetNumAtoms())
        ],
        "comment": comment,
        "num_bonds": xyz_mol.GetNumBonds(),
        "error": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/smiles.smi"),
        help="input SMILES file (default: tests/smiles.smi)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/xyz_read.jsonl"),
        help="output JSONL path (default: tests/golden/xyz_read.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = [build_record(smiles) for smiles in iter_smiles(args.input)]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
