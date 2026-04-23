#!/usr/bin/env python3
"""Generate RDKit V2000-body golden data for minimal molblock subset tests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from rdkit import Chem
from rdkit.Chem import AllChem


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        yield line


def molblock_body(block: str) -> str:
    lines = block.splitlines()
    if len(lines) < 4:
        return ""
    return "\n".join(lines[3:])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("tests/molblock_minimal_smiles.txt"),
        help="input smiles corpus (default: tests/molblock_minimal_smiles.txt)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("tests/golden/molblock_v2000_minimal.jsonl"),
        help="output jsonl (default: tests/golden/molblock_v2000_minimal.jsonl)",
    )
    args = parser.parse_args()

    rows = list(iter_smiles(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for smi in rows:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                rec = {
                    "smiles": smi,
                    "rdkit_ok": False,
                    "body": None,
                    "error": "MolFromSmiles returned None",
                }
            else:
                AllChem.Compute2DCoords(mol)
                mb = Chem.MolToMolBlock(mol)
                rec = {
                    "smiles": smi,
                    "rdkit_ok": True,
                    "body": molblock_body(mb),
                    "error": None,
                }
            f.write(json.dumps(rec, ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
