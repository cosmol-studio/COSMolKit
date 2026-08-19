#!/usr/bin/env python3
"""Generate exact Avalon bit records from the pinned RDKit adapter."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

from rdkit import Chem, RDLogger
from rdkit.Avalon import pyAvalonTools

PROFILE_PATH = Path(__file__).with_name("avalon_fingerprint_profile.json")
ENGINE_TAG = "AvalonToolkit_2.0.5-pre.3"
RDLogger.DisableLog("rdApp.*")


def iter_smiles(path: Path) -> Iterable[str]:
    for raw in path.read_text(encoding="utf-8").splitlines():
        value = raw.strip()
        if value and not value.startswith("#"):
            yield value


def load_profile() -> dict[str, object]:
    profile = json.loads(PROFILE_PATH.read_text(encoding="utf-8"))
    if profile.get("avalon_engine") != ENGINE_TAG:
        raise SystemExit("Avalon profile engine tag is not pinned")
    branches = profile.get("corpus_branches")
    focused = profile.get("focused_cases")
    if not isinstance(branches, list) or not branches:
        raise SystemExit("Avalon profile must contain corpus branches")
    if not isinstance(focused, list) or not focused:
        raise SystemExit("Avalon profile must contain focused cases")
    names = [branch.get("name") for branch in branches]
    if any(not isinstance(name, str) for name in names) or len(set(names)) != len(names):
        raise SystemExit("Avalon corpus branch names must be unique")
    return profile


def run_oracle(requests: list[tuple[str, dict[str, object]]]) -> list[dict[str, object]]:
    results = []
    for smiles, branch in requests:
        n_bits = int(branch["n_bits"])
        is_query = bool(branch["is_query"])
        bit_flags = int(str(branch["bit_flags"]), 16)
        try:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                raise ValueError("conversion")
            fingerprint = pyAvalonTools.GetAvalonFP(
                molecule, nBits=n_bits, isQuery=is_query, bitFlags=bit_flags
            )
            on_bits = [int(bit) for bit in fingerprint.GetOnBits()]
            bytes_value = bytearray(n_bits // 8)
            for bit in on_bits:
                if bit < len(bytes_value) * 8:
                    bytes_value[bit // 8] |= 1 << (bit % 8)
            results.append({
                "ok": True,
                "smiles": smiles,
                "n_bits": n_bits,
                "is_query": is_query,
                "bit_flags": f"0x{bit_flags:06x}",
                "bytes_hex": bytes(bytes_value).hex(),
                "counts": {},
                "on_bits": on_bits,
            })
        except Exception as error:
            results.append({"ok": False, "smiles": smiles, "error": str(error) or "conversion"})
    if len(results) != len(requests):
        raise SystemExit(f"Avalon oracle returned {len(results)} rows for {len(requests)} requests")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--focused-output", type=Path)
    args = parser.parse_args()
    profile = load_profile()
    branches = profile["corpus_branches"]
    focused = profile["focused_cases"]
    smiles_values = list(iter_smiles(args.input))
    corpus_requests = [(smiles, branch) for smiles in smiles_values for branch in branches]
    corpus_results = run_oracle(corpus_requests)
    rows = []
    offset = 0
    for row, smiles in enumerate(smiles_values):
        branch_records = {}
        for branch in branches:
            branch_records[branch["name"]] = {
                "parameters": branch,
                **corpus_results[offset],
            }
            offset += 1
        rows.append({"row": row, "smiles": smiles, "branches": branch_records})

    focused_requests = [(case["smiles"], case) for case in focused]
    focused_results = run_oracle(focused_requests)
    focused_rows = []
    for case, result in zip(focused, focused_results):
        focused_rows.append({"case": case, "result": result})

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")))
            handle.write("\n")
    if args.focused_output is not None:
        args.focused_output.parent.mkdir(parents=True, exist_ok=True)
        args.focused_output.write_text(
            json.dumps(
                {
                    "engine": ENGINE_TAG,
                    "profile": profile,
                    "cases": focused_rows,
                },
                sort_keys=True,
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )

    print(f"Wrote {len(rows)} corpus rows and {len(focused_rows)} focused cases")


if __name__ == "__main__":
    main()
