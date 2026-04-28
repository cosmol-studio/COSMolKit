#!/usr/bin/env python3
"""Generate RDKit ETKDG geometry golden data for tetrahedral stereo checks."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

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


def build_record(smiles: str) -> dict[str, object] | None:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return None

    chiral_centers = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if str(atom.GetChiralTag()) in {"CHI_TETRAHEDRAL_CW", "CHI_TETRAHEDRAL_CCW"}
    ]
    if not chiral_centers:
        return None

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    res = AllChem.EmbedMolecule(mol, params)
    if res != 0:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "centers": chiral_centers,
            "positions": None,
            "error": "EmbedMolecule failed",
        }

    conf = mol.GetConformer()
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "centers": chiral_centers,
        "positions": conf.GetPositions().tolist(),
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
        default=Path("tests/golden/tetrahedral_stereo_geometry.jsonl"),
        help="output JSONL path (default: tests/golden/tetrahedral_stereo_geometry.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = [record for smiles in iter_smiles(args.input) if (record := build_record(smiles))]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with args.output.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")

    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
