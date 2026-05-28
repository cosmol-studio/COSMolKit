#!/usr/bin/env python3
"""Generate RDKit XYZ reader golden data from RDKit built-in fixtures."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

from rdkit import Chem, RDLogger

EXPECTED_RDKIT_VERSION = "2026.3.1"
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "xyz" / "rdkit"
DEFAULT_OUTPUT = REPO_ROOT / "tests" / "golden" / "xyz_builtin_read.jsonl"

FIXTURE_FILES = [
    "acetate.xyz",
    "empty.xyz",
    "ethane.xyz",
    "nonexistant.xyz",
]


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="accepted for gen_all_rdkit_goldens.py compatibility; ignored",
    )
    parser.add_argument(
        "--fixture-dir",
        type=Path,
        default=DEFAULT_FIXTURE_DIR,
        help=f"directory containing copied RDKit XYZ fixtures (default: {DEFAULT_FIXTURE_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output JSONL path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def build_record(fixture_dir: Path, fixture: str) -> dict[str, object]:
    xyz_block = (fixture_dir / fixture).read_text(encoding="utf-8")
    try:
        mol = Chem.MolFromXYZBlock(xyz_block)
        if mol is None:
            return {
                "fixture": fixture,
                "rdkit_ok": False,
                "xyz_block": xyz_block,
                "atomic_numbers": None,
                "coordinates": None,
                "comment": None,
                "num_bonds": None,
                "error": "MolFromXYZBlock returned None",
            }
        conf = mol.GetConformer() if mol.GetNumConformers() else None
        return {
            "fixture": fixture,
            "rdkit_ok": True,
            "xyz_block": xyz_block,
            "atomic_numbers": [atom.GetAtomicNum() for atom in mol.GetAtoms()],
            "coordinates": [
                [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                for i in range(mol.GetNumAtoms())
            ]
            if conf is not None
            else [],
            "comment": mol.GetProp("_FileComments")
            if mol.HasProp("_FileComments")
            else None,
            "num_bonds": mol.GetNumBonds(),
            "error": None,
        }
    except Exception as err:
        return {
            "fixture": fixture,
            "rdkit_ok": False,
            "xyz_block": xyz_block,
            "atomic_numbers": None,
            "coordinates": None,
            "comment": None,
            "num_bonds": None,
            "error": str(err),
        }


def main() -> None:
    args = parse_args()
    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")

    records = [build_record(args.fixture_dir, fixture) for fixture in FIXTURE_FILES]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True, sort_keys=True) + "\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
