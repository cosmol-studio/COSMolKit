#!/usr/bin/env python3
"""Generate RDKit PrepareMolForDrawing golden data for tests/smiles.smi."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D

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
            "atoms": None,
            "bonds": None,
            "error": "MolFromSmiles returned None",
        }
    try:
        prepared = rdMolDraw2D.PrepareMolForDrawing(
            mol,
            kekulize=True,
            addChiralHs=True,
            wedgeBonds=True,
            forceCoords=True,
        )
        conf = prepared.GetConformer()
        atoms = []
        for atom in prepared.GetAtoms():
            pos = conf.GetAtomPosition(atom.GetIdx())
            atoms.append(
                {
                    "idx": atom.GetIdx(),
                    "atomic_num": atom.GetAtomicNum(),
                    "x": pos.x,
                    "y": pos.y,
                }
            )
        bonds = []
        for bond in prepared.GetBonds():
            bonds.append(
                {
                    "idx": bond.GetIdx(),
                    "begin": bond.GetBeginAtomIdx(),
                    "end": bond.GetEndAtomIdx(),
                    "bond_type": str(bond.GetBondType()),
                    "is_aromatic": bond.GetIsAromatic(),
                    "dir": str(bond.GetBondDir()),
                }
            )
        return {
            "smiles": smiles,
            "rdkit_ok": True,
            "atoms": atoms,
            "bonds": bonds,
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "atoms": None,
            "bonds": None,
            "error": f"{type(exc).__name__}: {exc}",
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
        default=Path("tests/golden/prepared_draw_molecule.jsonl"),
        help="output JSONL path (default: tests/golden/prepared_draw_molecule.jsonl)",
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
