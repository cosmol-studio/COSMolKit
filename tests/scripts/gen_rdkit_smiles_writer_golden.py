#!/usr/bin/env python3
"""Generate RDKit MolToSmiles golden data across writer parameter branches."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger

EXPECTED_RDKIT_VERSION = "2026.3.1"

BRANCHES = [
    {"name": "iso_true_kek_false_can_true", "isomericSmiles": True, "kekuleSmiles": False, "canonical": True},
    {"name": "iso_false_kek_false_can_true", "isomericSmiles": False, "kekuleSmiles": False, "canonical": True},
    {"name": "iso_true_kek_false_can_false", "isomericSmiles": True, "kekuleSmiles": False, "canonical": False},
    {"name": "iso_false_kek_false_can_false", "isomericSmiles": False, "kekuleSmiles": False, "canonical": False},
    {"name": "iso_true_kek_true_can_true", "isomericSmiles": True, "kekuleSmiles": True, "canonical": True},
    {"name": "iso_false_kek_true_can_true", "isomericSmiles": False, "kekuleSmiles": True, "canonical": True},
    {"name": "iso_true_kek_true_can_false", "isomericSmiles": True, "kekuleSmiles": True, "canonical": False},
    {"name": "iso_false_kek_true_can_false", "isomericSmiles": False, "kekuleSmiles": True, "canonical": False},
]


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


def branch_result(mol: Chem.Mol, branch: dict[str, object]) -> dict[str, object]:
    kwargs = {
        "isomericSmiles": branch["isomericSmiles"],
        "kekuleSmiles": branch["kekuleSmiles"],
        "canonical": branch["canonical"],
    }
    try:
        return {
            "ok": True,
            "smiles": Chem.MolToSmiles(mol, **kwargs),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "ok": False,
            "smiles": None,
            "error": f"{type(exc).__name__}: {exc}",
        }


def build_record(smiles: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "branches": {},
            "error": "MolFromSmiles returned None",
        }
    branches = {branch["name"]: branch_result(mol, branch) for branch in BRANCHES}
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "branches": branches,
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
        default=Path("tests/golden/smiles_writer.jsonl"),
        help="output JSONL path (default: tests/golden/smiles_writer.jsonl)",
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
