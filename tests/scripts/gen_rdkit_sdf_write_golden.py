#!/usr/bin/env python3
"""Generate RDKit MolToMolBlock writer golden data across parameter branches."""

from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing
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


def render_source_molblock(source: Chem.Mol | None) -> tuple[str | None, str | None]:
    if source is None:
        return None, "source molecule generation failed"
    try:
        return (
            Chem.MolToMolBlock(
                source, includeStereo=True, kekulize=True, forceV3000=True
            ),
            None,
        )
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def parse_source_molblock(block: str | None) -> tuple[Chem.Mol | None, str | None]:
    if block is None:
        return None, "source molblock generation failed"
    try:
        parsed = Chem.MolFromMolBlock(block, sanitize=True, removeHs=False)
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"
    if parsed is None:
        return None, "MolFromMolBlock returned None"
    return parsed, None


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
    source_2d_molblock, source_2d_molblock_error = render_source_molblock(source_2d)
    source_3d_molblock, source_3d_molblock_error = render_source_molblock(source_3d)
    source_2d_for_write, source_2d_parse_error = parse_source_molblock(source_2d_molblock)
    source_3d_for_write, source_3d_parse_error = parse_source_molblock(source_3d_molblock)
    source_2d_error = error_2d or source_2d_molblock_error or source_2d_parse_error
    source_3d_error = error_3d or source_3d_molblock_error or source_3d_parse_error
    branches = {}
    for dimension, source, source_error in [
        ("2d", source_2d_for_write, source_2d_error),
        ("3d", source_3d_for_write, source_3d_error),
    ]:
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
        "source_2d_error": source_2d_error,
        "source_3d_molblock": source_3d_molblock,
        "source_3d_error": source_3d_error,
        "branches": branches,
        "error": None,
    }


def build_indexed_record(item: tuple[int, str]) -> tuple[int, dict[str, object]]:
    row_idx, smiles = item
    return row_idx, build_record(smiles)


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
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="parallel RDKit worker processes (default: 1)",
    )
    args = parser.parse_args()

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    smiles_list = list(iter_smiles(args.input))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp_output = args.output.with_suffix(args.output.suffix + ".tmp")
    items = list(enumerate(smiles_list, start=1))

    with tmp_output.open("w", encoding="utf-8") as handle:
        next_to_write = 1
        pending: dict[int, dict[str, object]] = {}
        completed = 0
        if args.jobs <= 1:
            iterator = map(build_indexed_record, items)
            for row_idx, record in iterator:
                completed += 1
                if completed == 1 or completed % 50 == 0 or completed == len(smiles_list):
                    print(f"[sdf-write] completed {completed}/{len(smiles_list)}", flush=True)
                pending[row_idx] = record
                while next_to_write in pending:
                    record_to_write = pending.pop(next_to_write)
                    handle.write(json.dumps(record_to_write, ensure_ascii=True))
                    handle.write("\n")
                    next_to_write += 1
                handle.flush()
        else:
            with multiprocessing.Pool(processes=args.jobs) as pool:
                iterator = pool.imap_unordered(build_indexed_record, items, chunksize=1)
                for row_idx, record in iterator:
                    completed += 1
                    if completed == 1 or completed % 50 == 0 or completed == len(smiles_list):
                        print(
                            f"[sdf-write] completed {completed}/{len(smiles_list)}",
                            flush=True,
                        )
                    pending[row_idx] = record
                    while next_to_write in pending:
                        record_to_write = pending.pop(next_to_write)
                        handle.write(json.dumps(record_to_write, ensure_ascii=True))
                        handle.write("\n")
                        next_to_write += 1
                    handle.flush()
    tmp_output.replace(args.output)
    print(f"Wrote {len(smiles_list)} records to {args.output}")


if __name__ == "__main__":
    main()
