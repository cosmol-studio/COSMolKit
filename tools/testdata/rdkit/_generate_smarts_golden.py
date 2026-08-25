#!/usr/bin/env python3
"""Generate pinned-RDKit SMARTS parse, write, and match golden data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def build_record(case: dict[str, Any]) -> dict[str, Any]:
    smarts = str(case["smarts"])
    query = Chem.MolFromSmarts(smarts)
    record: dict[str, Any] = {
        "case_id": str(case["case_id"]),
        "kind": str(case["kind"]),
        "source": str(case["source"]),
        "smarts": smarts,
        "target_smiles": case.get("target_smiles"),
        "parse_ok": query is not None,
        "num_atoms": None,
        "num_bonds": None,
        "written_smarts": None,
        "atom_mappings": None,
    }
    if query is None:
        return record

    record["num_atoms"] = query.GetNumAtoms()
    record["num_bonds"] = query.GetNumBonds()
    record["written_smarts"] = Chem.MolToSmarts(query)

    target_smiles = case.get("target_smiles")
    if target_smiles is not None:
        target = Chem.MolFromSmiles(str(target_smiles))
        if target is None:
            raise ValueError(f"{case['case_id']}: RDKit rejected target SMILES")
        record["atom_mappings"] = [
            list(mapping) for mapping in target.GetSubstructMatches(query)
        ]
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    document = json.loads(args.input.read_text(encoding="utf-8"))
    cases = document.get("cases")
    if not isinstance(cases, list):
        raise ValueError("SMARTS corpus must contain a top-level cases array")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for case in cases:
            if not isinstance(case, dict):
                raise ValueError("every SMARTS corpus case must be an object")
            handle.write(json.dumps(build_record(case), ensure_ascii=True))
            handle.write("\n")
    print(f"Wrote {len(cases)} SMARTS records to {args.output}")


if __name__ == "__main__":
    main()
