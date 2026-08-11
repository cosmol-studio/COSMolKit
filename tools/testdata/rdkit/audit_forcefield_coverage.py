#!/usr/bin/env python3
"""Prepare a bounded ChEMBL force-field coverage audit corpus and RDKit oracle."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger

from _generate_forcefield_coverage_golden import build_record
from _generate_forcefield_params_golden import assert_rdkit_version


def read_jsonl(path: Path):
    with path.open(encoding="utf-8") as source:
        for line_number, raw in enumerate(source, start=1):
            if not raw.strip():
                continue
            record: Any = json.loads(raw)
            if not isinstance(record, dict) or not isinstance(record.get("smiles"), str):
                raise ValueError(f"{path}:{line_number}: expected an object with string smiles")
            yield record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-jsonl", type=Path, required=True)
    parser.add_argument("--smiles-output", type=Path, required=True)
    parser.add_argument("--golden-output", type=Path, required=True)
    parser.add_argument("--max-atoms", type=int, default=80)
    parser.add_argument("--scan-limit", type=int, default=None)
    args = parser.parse_args()

    if args.max_atoms < 1:
        parser.error("--max-atoms must be >= 1")

    assert_rdkit_version()
    RDLogger.DisableLog("rdApp.*")
    args.smiles_output.parent.mkdir(parents=True, exist_ok=True)
    args.golden_output.parent.mkdir(parents=True, exist_ok=True)

    scanned = 0
    selected = 0
    parse_failures = 0
    oversize = 0
    with args.smiles_output.open("w", encoding="utf-8") as smiles_handle, args.golden_output.open(
        "w", encoding="utf-8"
    ) as golden_handle:
        for source_record in read_jsonl(args.input_jsonl):
            if args.scan_limit is not None and scanned >= args.scan_limit:
                break
            scanned += 1
            smiles = source_record["smiles"]
            mol = Chem.MolFromSmiles(smiles, sanitize=True)
            if mol is None:
                parse_failures += 1
                continue
            if mol.GetNumAtoms() > args.max_atoms:
                oversize += 1
                continue

            smiles_handle.write(smiles)
            smiles_handle.write("\n")
            golden_handle.write(json.dumps(build_record(smiles), ensure_ascii=True))
            golden_handle.write("\n")
            selected += 1

    print(
        json.dumps(
            {
                "input": str(args.input_jsonl),
                "max_atoms": args.max_atoms,
                "parse_failures": parse_failures,
                "scanned": scanned,
                "selected": selected,
                "skipped_oversize": oversize,
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
