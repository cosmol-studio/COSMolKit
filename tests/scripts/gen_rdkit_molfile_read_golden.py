#!/usr/bin/env python3
"""Generate RDKit molfile read-parity golden data.

This uses Chem.MolFromMolBlock() directly, not SDMolSupplier, so the golden
documents single-record .mol behavior without SDF record separators or data
fields. The case matrix mirrors SDF read parity: 2D/3D coordinates, V2000/V3000,
and explicit stereo markers versus coordinate-only records.
"""

from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger

from gen_rdkit_sdf_read_golden import (
    EXPECTED_RDKIT_VERSION,
    atom_features,
    bond_features,
    chiral_tags,
    iter_smiles,
    make_2d_source,
    make_3d_source,
    positions,
    render_molblock,
    smiles_pair,
    strip_stereo_markers,
)


def assert_rdkit_version() -> None:
    actual = version("rdkit")
    assert actual == EXPECTED_RDKIT_VERSION, (
        f"RDKit version mismatch: expected {EXPECTED_RDKIT_VERSION}, got {actual}"
    )


def parse_molblock(block: str) -> tuple[Chem.Mol | None, str | None]:
    try:
        parsed = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
    except Exception as exc:
        return None, str(exc)
    if parsed is None:
        return None, "MolFromMolBlock returned None"
    return parsed, None


def build_case(
    smiles: str,
    source: Chem.Mol | None,
    source_error: str | None,
    *,
    dimension: str,
    fmt: str,
    markers: str,
) -> dict[str, object]:
    case_id = f"{dimension.lower()}_{fmt.lower()}_{markers}"
    base = {
        "smiles": smiles,
        "case_id": case_id,
        "dimension": dimension,
        "format": fmt,
        "stereo_markers": markers,
    }
    if source is None:
        return {
            **base,
            "rdkit_ok": False,
            "molblock": None,
            "error": source_error or "source molecule generation failed",
        }

    block, render_error = render_molblock(source, fmt)
    if block is None:
        return {
            **base,
            "rdkit_ok": False,
            "molblock": None,
            "error": render_error,
        }
    if markers == "coords_only":
        block = strip_stereo_markers(block, fmt)

    parsed, parse_error = parse_molblock(block)
    if parsed is None:
        return {
            **base,
            "rdkit_ok": False,
            "molblock": block,
            "error": parse_error,
        }

    return {
        **base,
        "rdkit_ok": True,
        "molblock": block,
        "atoms": atom_features(parsed),
        "bonds": bond_features(parsed),
        "positions": positions(parsed),
        "chiral_tags": chiral_tags(parsed),
        "smiles_out": smiles_pair(parsed),
        "error": None,
    }


def build_records(smiles: str) -> list[dict[str, object]]:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    if mol is None:
        return [
            {
                "smiles": smiles,
                "case_id": f"{dimension.lower()}_{fmt.lower()}_{markers}",
                "dimension": dimension,
                "format": fmt,
                "stereo_markers": markers,
                "rdkit_ok": False,
                "molblock": None,
                "error": "MolFromSmiles returned None",
            }
            for dimension in ["2D", "3D"]
            for fmt in ["V2000", "V3000"]
            for markers in ["with_markers", "coords_only"]
        ]

    source_2d, error_2d = make_2d_source(mol)
    source_3d, error_3d = make_3d_source(mol)
    records: list[dict[str, object]] = []
    for dimension, source, source_error in [
        ("2D", source_2d, error_2d),
        ("3D", source_3d, error_3d),
    ]:
        for fmt in ["V2000", "V3000"]:
            for markers in ["with_markers", "coords_only"]:
                records.append(
                    build_case(
                        smiles,
                        source,
                        source_error,
                        dimension=dimension,
                        fmt=fmt,
                        markers=markers,
                    )
                )
    return records


def build_all_records(smiles_rows: Iterable[str]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for smiles in smiles_rows:
        rows.extend(build_records(smiles))
    return rows


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
        default=Path("tests/golden/molfile_read.jsonl"),
        help="output JSONL path (default: tests/golden/molfile_read.jsonl)",
    )
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    smiles_rows = list(iter_smiles(args.input))
    if args.limit is not None:
        smiles_rows = smiles_rows[: args.limit]
    rows = build_all_records(smiles_rows)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for record in rows:
            handle.write(json.dumps(record, ensure_ascii=True))
            handle.write("\n")
    print(f"Wrote {len(rows)} records to {args.output}")


if __name__ == "__main__":
    main()
