#!/usr/bin/env python3
"""Generate RDKit MolToMolBlock writer golden data across parameter branches."""

from __future__ import annotations

import argparse
import itertools
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


def make_2d_source(mol: Chem.Mol) -> tuple[Chem.Mol | None, str | None]:
    source = Chem.Mol(mol)
    try:
        AllChem.Compute2DCoords(source)
        if source.GetNumConformers():
            Chem.WedgeMolBonds(source, source.GetConformer())
        return source, None
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def make_3d_source(mol: Chem.Mol) -> tuple[Chem.Mol | None, str | None]:
    source = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    try:
        status = AllChem.EmbedMolecule(source, params)
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"
    if status != 0:
        return None, f"EmbedMolecule failed with status {status}"
    return source, None


def molblock_body(block: str) -> str:
    lines = block.splitlines()
    if len(lines) < 4:
        return ""
    return "\n".join(lines[3:])


def branch_result(
    source: Chem.Mol | None,
    source_error: str | None,
    *,
    dimension: str,
    include_stereo: bool,
    kekulize: bool,
    force_v3000: bool,
) -> dict[str, object]:
    params = {
        "dimension": dimension,
        "include_stereo": include_stereo,
        "kekulize": kekulize,
        "force_v3000": force_v3000,
    }
    if source is None:
        return {
            "params": params,
            "ok": False,
            "body": None,
            "error": source_error or "source molecule generation failed",
        }
    try:
        block = Chem.MolToMolBlock(
            source,
            includeStereo=include_stereo,
            confId=-1,
            kekulize=kekulize,
            forceV3000=force_v3000,
        )
        return {
            "params": params,
            "ok": True,
            "body": molblock_body(block),
            "error": None,
        }
    except Exception as exc:  # noqa: BLE001
        return {
            "params": params,
            "ok": False,
            "body": None,
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
    source_2d, error_2d = make_2d_source(mol)
    source_3d, error_3d = make_3d_source(mol)
    source_2d_molblock = (
        Chem.MolToMolBlock(source_2d, includeStereo=True, kekulize=True, forceV3000=True)
        if source_2d is not None
        else None
    )
    source_3d_molblock = (
        Chem.MolToMolBlock(source_3d, includeStereo=True, kekulize=True, forceV3000=True)
        if source_3d is not None
        else None
    )
    source_2d_for_write = (
        Chem.MolFromMolBlock(source_2d_molblock, sanitize=True, removeHs=False)
        if source_2d_molblock is not None
        else None
    )
    source_3d_for_write = (
        Chem.MolFromMolBlock(source_3d_molblock, sanitize=True, removeHs=False)
        if source_3d_molblock is not None
        else None
    )
    branches = {}
    for dimension, source, source_error in [
        ("2d", source_2d, error_2d),
        ("3d", source_3d, error_3d),
    ]:
        source = source_2d_for_write if dimension == "2d" else source_3d_for_write
        for include_stereo, kekulize, force_v3000 in itertools.product(
            [False, True], repeat=3
        ):
            name = (
                f"{dimension}_stereo{int(include_stereo)}_"
                f"kek{int(kekulize)}_v3k{int(force_v3000)}"
            )
            branches[name] = branch_result(
                source,
                source_error,
                dimension=dimension,
                include_stereo=include_stereo,
                kekulize=kekulize,
                force_v3000=force_v3000,
            )
    return {
        "smiles": smiles,
        "rdkit_ok": True,
        "source_2d_molblock": source_2d_molblock,
        "source_3d_molblock": source_3d_molblock,
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
        default=Path("tests/golden/sdf_write.jsonl"),
        help="output JSONL path (default: tests/golden/sdf_write.jsonl)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    records = [build_record(smiles) for smiles in iter_smiles(args.input)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in records:
            handle.write(json.dumps(record, ensure_ascii=True))
            handle.write("\n")
    print(f"Wrote {len(records)} records to {args.output}")


if __name__ == "__main__":
    main()
