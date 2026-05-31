#!/usr/bin/env python3
"""Generate RDKit single-conformer library golden data from tests/smiles.smi."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDistGeom

EXPECTED_RDKIT_VERSION = "2026.3.1"
CONFORMER_LIBRARY_SEED = 61453
# Keep library parity on a deterministic embed-attempt budget. RDKit's
# timeout is wall-clock based, so it is unsuitable for CI-stable parity rows.
CONFORMER_LIBRARY_MAX_ITERATIONS = 3
CONFORMER_LIBRARY_TIMEOUT = 0


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
            "seed": CONFORMER_LIBRARY_SEED,
            "preset": "ETKDGv3",
            "max_iterations": CONFORMER_LIBRARY_MAX_ITERATIONS,
            "timeout": CONFORMER_LIBRARY_TIMEOUT,
            "rdkit_parse_ok": False,
            "rdkit_add_hs_ok": False,
            "rdkit_embed_ok": False,
            "status": None,
            "coords": None,
            "error_stage": "parse",
            "error": "MolFromSmiles returned None",
        }

    try:
        mol = Chem.AddHs(mol)
    except Exception as err:
        return {
            "smiles": smiles,
            "seed": CONFORMER_LIBRARY_SEED,
            "preset": "ETKDGv3",
            "max_iterations": CONFORMER_LIBRARY_MAX_ITERATIONS,
            "timeout": CONFORMER_LIBRARY_TIMEOUT,
            "rdkit_parse_ok": True,
            "rdkit_add_hs_ok": False,
            "rdkit_embed_ok": False,
            "status": None,
            "coords": None,
            "error_stage": "add_hs",
            "error": str(err),
        }

    try:
        params = rdDistGeom.ETKDGv3()
        params.maxIterations = CONFORMER_LIBRARY_MAX_ITERATIONS
        params.randomSeed = CONFORMER_LIBRARY_SEED
        params.numThreads = 1
        params.timeout = CONFORMER_LIBRARY_TIMEOUT
        status = int(rdDistGeom.EmbedMolecule(mol, params))
        if status != 0:
            return {
                "smiles": smiles,
                "seed": CONFORMER_LIBRARY_SEED,
                "preset": "ETKDGv3",
                "max_iterations": CONFORMER_LIBRARY_MAX_ITERATIONS,
                "timeout": CONFORMER_LIBRARY_TIMEOUT,
                "rdkit_parse_ok": True,
                "rdkit_add_hs_ok": True,
                "rdkit_embed_ok": False,
                "status": status,
                "coords": None,
                "error_stage": "embed",
                "error": f"EmbedMolecule returned {status}",
            }
        conf = mol.GetConformer(status)
        coords = [
            [float(conf.GetAtomPosition(i)[axis]) for axis in range(3)]
            for i in range(mol.GetNumAtoms())
        ]
        return {
            "smiles": smiles,
            "seed": CONFORMER_LIBRARY_SEED,
            "preset": "ETKDGv3",
            "max_iterations": CONFORMER_LIBRARY_MAX_ITERATIONS,
            "timeout": CONFORMER_LIBRARY_TIMEOUT,
            "rdkit_parse_ok": True,
            "rdkit_add_hs_ok": True,
            "rdkit_embed_ok": True,
            "status": status,
            "coords": coords,
            "error_stage": None,
            "error": None,
        }
    except Exception as err:
        return {
            "smiles": smiles,
            "seed": CONFORMER_LIBRARY_SEED,
            "preset": "ETKDGv3",
            "max_iterations": CONFORMER_LIBRARY_MAX_ITERATIONS,
            "timeout": CONFORMER_LIBRARY_TIMEOUT,
            "rdkit_parse_ok": True,
            "rdkit_add_hs_ok": True,
            "rdkit_embed_ok": False,
            "status": None,
            "coords": None,
            "error_stage": "embed",
            "error": str(err),
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
        default=Path("tests/golden/conformer_generation_library.jsonl"),
        help="output JSONL path (default: tests/golden/conformer_generation_library.jsonl)",
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
