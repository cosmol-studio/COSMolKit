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
import tempfile
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


def parse_molblock_with_params(
    block: str, *, sanitize: bool, remove_hs: bool, strict_parsing: bool
) -> tuple[Chem.Mol | None, str | None]:
    try:
        parsed = Chem.MolFromMolBlock(
            block,
            sanitize=sanitize,
            removeHs=remove_hs,
            strictParsing=strict_parsing,
        )
    except Exception as exc:
        return None, str(exc)
    if parsed is None:
        return None, "MolFromMolBlock returned None"
    return parsed, None


def parse_molfile_with_params(
    block: str, *, sanitize: bool, remove_hs: bool, strict_parsing: bool
) -> tuple[Chem.Mol | None, str | None]:
    with tempfile.NamedTemporaryFile("w", suffix=".mol", encoding="utf-8") as handle:
        handle.write(block)
        handle.flush()
        try:
            parsed = Chem.MolFromMolFile(
                handle.name,
                sanitize=sanitize,
                removeHs=remove_hs,
                strictParsing=strict_parsing,
            )
        except Exception as exc:
            return None, str(exc)
    if parsed is None:
        return None, "MolFromMolFile returned None"
    return parsed, None


def delayed_sanitize(mol: Chem.Mol | None) -> tuple[Chem.Mol | None, str | None]:
    if mol is None:
        return None, "source molecule is None"
    copy = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(copy)
    except Exception as exc:
        return None, str(exc)
    return copy, None


def delayed_remove_hs(mol: Chem.Mol | None) -> tuple[Chem.Mol | None, str | None]:
    if mol is None:
        return None, "source molecule is None"
    try:
        removed = Chem.RemoveHs(Chem.Mol(mol))
    except Exception as exc:
        return None, str(exc)
    if removed is None:
        return None, "RemoveHs returned None"
    return removed, None


def mol_summary(mol: Chem.Mol) -> dict[str, object]:
    return {
        "atoms": atom_features(mol),
        "bonds": bond_features(mol),
        "positions": positions(mol),
        "chiral_tags": chiral_tags(mol),
        "smiles_out": smiles_pair(mol),
    }


def error_row(base: dict[str, object], block: str | None, error: str | None) -> dict[str, object]:
    return {
        **base,
        "rdkit_ok": False,
        "molblock": block,
        "error": error or "RDKit operation failed",
    }


def success_row(base: dict[str, object], block: str, mol: Chem.Mol) -> dict[str, object]:
    return {
        **base,
        "rdkit_ok": True,
        "molblock": block,
        **mol_summary(mol),
        "error": None,
    }


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
        "api": "MolFromMolBlock",
        "operation": "read",
        "sanitize": True,
        "remove_hs": False,
        "strict_parsing": True,
    }
    if source is None:
        return error_row(base, None, source_error or "source molecule generation failed")

    block, render_error = render_molblock(source, fmt)
    if block is None:
        return error_row(base, None, render_error)
    if markers == "coords_only":
        block = strip_stereo_markers(block, fmt)

    parsed, parse_error = parse_molblock(block)
    if parsed is None:
        return error_row(base, block, parse_error)

    return success_row(base, block, parsed)


def build_parameterized_cases(
    smiles: str,
    block: str,
    *,
    dimension: str,
    fmt: str,
    markers: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for api_name, parser_fn in [
        ("MolFromMolBlock", parse_molblock_with_params),
        ("MolFromMolFile", parse_molfile_with_params),
    ]:
        for sanitize in [True, False]:
            for remove_hs in [True, False]:
                for strict_parsing in [True, False]:
                    case_id = (
                        f"{dimension.lower()}_{fmt.lower()}_{markers}_"
                        f"{api_name.lower()}_sanitize_{int(sanitize)}_"
                        f"removehs_{int(remove_hs)}_strict_{int(strict_parsing)}"
                    )
                    base = {
                        "smiles": smiles,
                        "case_id": case_id,
                        "dimension": dimension,
                        "format": fmt,
                        "stereo_markers": markers,
                        "api": api_name,
                        "operation": "read",
                        "sanitize": sanitize,
                        "remove_hs": remove_hs,
                        "strict_parsing": strict_parsing,
                    }
                    parsed, parse_error = parser_fn(
                        block,
                        sanitize=sanitize,
                        remove_hs=remove_hs,
                        strict_parsing=strict_parsing,
                    )
                    if parsed is None:
                        rows.append(error_row(base, block, parse_error))
                        continue
                    rows.append(success_row(base, block, parsed))

                    if not sanitize:
                        delayed, delayed_error = delayed_sanitize(parsed)
                        delayed_base = {
                            **base,
                            "case_id": f"{case_id}_delayed_sanitize",
                            "operation": "delayed_sanitize",
                        }
                        rows.append(
                            error_row(delayed_base, block, delayed_error)
                            if delayed is None
                            else success_row(delayed_base, block, delayed)
                        )
                    if not remove_hs:
                        delayed, delayed_error = delayed_remove_hs(parsed)
                        delayed_base = {
                            **base,
                            "case_id": f"{case_id}_delayed_remove_hs",
                            "operation": "delayed_remove_hs",
                        }
                        rows.append(
                            error_row(delayed_base, block, delayed_error)
                            if delayed is None
                            else success_row(delayed_base, block, delayed)
                        )
    return rows


def build_failure_cases() -> list[dict[str, object]]:
    cases = [
        ("empty_input", ""),
        ("missing_counts_line", "missing counts\n  COSMolKit      2D\ncomment\n"),
        (
            "short_counts_line",
            "short counts\n  COSMolKit      2D\ncomment\n  1\nM  END\n",
        ),
        (
            "invalid_counts_line",
            "bad counts\n  COSMolKit      2D\ncomment\nxxx  0  0  0  0  0  0  0  0  0999 V2000\nM  END\n",
        ),
        (
            "missing_m_end",
            "missing end\n  COSMolKit      2D\ncomment\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n",
        ),
    ]
    rows: list[dict[str, object]] = []
    for label, block in cases:
        for api_name, parser_fn in [
            ("MolFromMolBlock", parse_molblock_with_params),
            ("MolFromMolFile", parse_molfile_with_params),
        ]:
            for strict_parsing in [True, False]:
                base = {
                    "smiles": None,
                    "case_id": f"failure_{label}_{api_name.lower()}_strict_{int(strict_parsing)}",
                    "dimension": "invalid",
                    "format": "invalid",
                    "stereo_markers": "invalid",
                    "api": api_name,
                    "operation": "failure",
                    "sanitize": True,
                    "remove_hs": True,
                    "strict_parsing": strict_parsing,
                }
                parsed, parse_error = parser_fn(
                    block,
                    sanitize=True,
                    remove_hs=True,
                    strict_parsing=strict_parsing,
                )
                rows.append(
                    error_row(base, block, parse_error)
                    if parsed is None
                    else success_row(base, block, parsed)
                )
    return rows


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
                legacy = build_case(
                    smiles,
                    source,
                    source_error,
                    dimension=dimension,
                    fmt=fmt,
                    markers=markers,
                )
                records.append(legacy)
                block = legacy.get("molblock")
                if isinstance(block, str):
                    records.extend(
                        build_parameterized_cases(
                            smiles,
                            block,
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
    rows.extend(build_failure_cases())
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
