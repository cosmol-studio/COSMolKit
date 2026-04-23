#!/usr/bin/env python3
"""Generate RDKit atomic-number golden data from tests/smiles.txt."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from rdkit import Chem


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
            "atomic_nums": None,
            "error": "MolFromSmiles returned None",
        }

    atomic_nums = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "atomic_nums": atomic_nums,
        "error": None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/smiles.txt"),
        help="input SMILES file (default: tests/smiles.txt)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/atomic_nums.jsonl"),
        help="output JSONL path (default: tests/golden/atomic_nums.jsonl)",
    )
    args = parser.parse_args()

    smiles_rows = list(iter_smiles(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for smiles in smiles_rows:
            record = build_record(smiles)
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(smiles_rows)} records to {args.output}")


if __name__ == "__main__":
    main()

