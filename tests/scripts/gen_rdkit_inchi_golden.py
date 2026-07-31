#!/usr/bin/env python3
"""Generate pinned-RDKit MolToInchi golden data for a SMILES corpus."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

from rdkit import Chem, RDLogger


EXPECTED_RDKIT_VERSION = "2026.3.1"


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def load_smiles(path: Path) -> list[str]:
    return [
        line
        for raw in path.read_text(encoding="utf-8").splitlines()
        if (line := raw.strip()) and not line.startswith("#")
    ]


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
        default=Path("tests/golden/smiles_small/inchi.jsonl"),
        help="output JSONL path",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    smiles_rows = load_smiles(args.input)
    records: list[dict[str, object]] = []

    for row, smiles in enumerate(smiles_rows, start=1):
        molecule = Chem.MolFromSmiles(smiles, sanitize=True)
        if molecule is None:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "inchi": None,
                    "error": "MolFromSmiles returned None",
                }
            )
            continue

        try:
            inchi = Chem.MolToInchi(molecule)
        except Exception as error:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "inchi": None,
                    "error": f"{type(error).__name__}: {error}",
                }
            )
            continue

        if not inchi:
            records.append(
                {
                    "row": row,
                    "smiles": smiles,
                    "rdkit_ok": False,
                    "inchi": None,
                    "error": "MolToInchi returned an empty string",
                }
            )
            continue

        records.append(
            {
                "row": row,
                "smiles": smiles,
                "rdkit_ok": True,
                "inchi": inchi,
                "error": None,
            }
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as output:
        for record in records:
            output.write(json.dumps(record, ensure_ascii=True))
            output.write("\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
