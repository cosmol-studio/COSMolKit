#!/usr/bin/env python3
"""Generate RDKit MMFF built-in fixture golden data."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

EXPECTED_RDKIT_VERSION = "2026.3.1"
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "forcefield" / "mmff" / "rdkit"
DEFAULT_OUTPUT = REPO_ROOT / "tests" / "golden" / "mmff_builtin.jsonl"

FIXTURE_FILES = [
    "MMFF94_dative.smi",
    "MMFF94_hypervalent.smi",
    "MMFF94s_dative.smi",
    "MMFF94s_hypervalent.smi",
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
        help=f"directory containing copied RDKit MMFF .smi fixtures (default: {DEFAULT_FIXTURE_DIR})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"output JSONL path (default: {DEFAULT_OUTPUT})",
    )
    return parser.parse_args()


def variant_for_fixture(fixture_name: str) -> str:
    if fixture_name.startswith("MMFF94s_"):
        return "MMFF94s"
    if fixture_name.startswith("MMFF94_"):
        return "MMFF94"
    raise ValueError(f"cannot infer MMFF variant from fixture name: {fixture_name}")


def iter_fixture_rows(path: Path) -> list[tuple[int, str, str]]:
    rows: list[tuple[int, str, str]] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line_number == 1 and line == "SMILES Name":
            continue
        parts = line.split()
        if len(parts) != 2:
            raise ValueError(
                f"{path} line {line_number}: expected RDKit 'SMILES Name' row, got {line!r}"
            )
        rows.append((line_number, parts[0], parts[1]))
    return rows


def mmff_result(smiles: str, variant: str) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "rdkit_ok": False,
            "num_atoms": None,
            "has_all": None,
            "props_ok": False,
            "atom_types": None,
            "error": "MolFromSmiles returned None",
        }
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=variant)
        atom_types = (
            [int(props.GetMMFFAtomType(i)) for i in range(mol.GetNumAtoms())]
            if props is not None
            else None
        )
        return {
            "rdkit_ok": True,
            "num_atoms": int(mol.GetNumAtoms()),
            "has_all": bool(AllChem.MMFFHasAllMoleculeParams(mol)),
            "props_ok": props is not None,
            "atom_types": atom_types,
            "error": None,
        }
    except Exception as err:
        return {
            "rdkit_ok": True,
            "num_atoms": int(mol.GetNumAtoms()),
            "has_all": None,
            "props_ok": False,
            "atom_types": None,
            "error": str(err),
        }


def build_records(fixture_dir: Path) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for fixture_name in FIXTURE_FILES:
        path = fixture_dir / fixture_name
        variant = variant_for_fixture(fixture_name)
        for line_number, smiles, row_name in iter_fixture_rows(path):
            record = {
                "fixture": fixture_name,
                "line_number": line_number,
                "row_name": row_name,
                "smiles": smiles,
                "variant": variant,
            }
            record.update(mmff_result(smiles, variant))
            records.append(record)
    return records


def main() -> None:
    args = parse_args()
    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")

    records = build_records(args.fixture_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
