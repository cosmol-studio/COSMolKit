#!/usr/bin/env python3
"""Generate RDKit Kekulize(clearAromaticFlags=False) golden data."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem

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


def bond_type_name(bond: Chem.Bond) -> str:
    return str(bond.GetBondType())


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "parse_ok": False,
            "parse_error": "MolFromSmiles returned None",
            "kekulize_ok": False,
            "kekulize_error": "parse failed",
            "atom_is_aromatic": None,
            "bond_types": None,
            "bond_is_aromatic": None,
        }

    try:
        Chem.Kekulize(mol, clearAromaticFlags=False)
    except Exception as exc:  # noqa: BLE001
        return {
            "smiles": smiles,
            "parse_ok": True,
            "parse_error": None,
            "kekulize_ok": False,
            "kekulize_error": f"{type(exc).__name__}: {exc}",
            "atom_is_aromatic": None,
            "bond_types": None,
            "bond_is_aromatic": None,
        }

    return {
        "smiles": smiles,
        "parse_ok": True,
        "parse_error": None,
        "kekulize_ok": True,
        "kekulize_error": None,
        "atom_is_aromatic": [atom.GetIsAromatic() for atom in mol.GetAtoms()],
        "bond_types": [bond_type_name(bond) for bond in mol.GetBonds()],
        "bond_is_aromatic": [bond.GetIsAromatic() for bond in mol.GetBonds()],
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
        default=Path("tests/golden/kekulize_clear_flags_false.jsonl"),
        help="output JSONL path (default: tests/golden/kekulize_clear_flags_false.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    rows = [build_record(smiles) for smiles in iter_smiles(args.input)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        for row in rows:
            f.write(json.dumps(row, ensure_ascii=True))
            f.write("\n")
    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
