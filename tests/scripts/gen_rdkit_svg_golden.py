#!/usr/bin/env python3
"""Generate RDKit MolDraw2DSVG golden data for tests/smiles.smi."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import rdMolDraw2D

EXPECTED_RDKIT_VERSION = "2026.3.1"
DEFAULT_WIDTH = 300
DEFAULT_HEIGHT = 300


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


def svg_for_mol(mol: Chem.Mol, width: int, height: int) -> str:
    # The final SVG text is the Python-exposed RDKit drawing boundary. The
    # noFreetype=True constructor path uses DrawTextSVG's built-in metrics,
    # which is the stable source-level target for the Rust implementation.
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height, -1, -1, True)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def build_record(smiles: str, width: int, height: int) -> dict[str, object]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "width": width,
            "height": height,
            "svg": None,
            "error": "MolFromSmiles returned None",
        }
    try:
        return {
            "smiles": smiles,
            "rdkit_ok": True,
            "width": width,
            "height": height,
            "svg": svg_for_mol(mol, width, height),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "smiles": smiles,
            "rdkit_ok": False,
            "width": width,
            "height": height,
            "svg": None,
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
        default=Path("tests/golden/svg_drawer.jsonl"),
        help="output JSONL path (default: tests/golden/svg_drawer.jsonl)",
    )
    parser.add_argument("--width", type=int, default=DEFAULT_WIDTH)
    parser.add_argument("--height", type=int, default=DEFAULT_HEIGHT)
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = [build_record(smiles, args.width, args.height) for smiles in iter_smiles(args.input)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as f:
        for record in records:
            f.write(json.dumps(record, ensure_ascii=True))
            f.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
