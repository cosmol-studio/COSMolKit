#!/usr/bin/env python3
"""Exact ChEMBL tautomer parity against the pinned RDKit enumerator."""

from __future__ import annotations

import argparse
import itertools
import json
import os
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Any

ORACLE_DIR = Path(__file__).resolve().parents[3] / "tools/testdata/rdkit"
if str(ORACLE_DIR) not in sys.path:
    sys.path.insert(0, str(ORACLE_DIR))
import _tautomer_oracle as oracle  # noqa: E402


class Audit:
    def __init__(self, output: Path, max_examples: int) -> None:
        self.output = output
        self.findings_path = output.with_suffix(".findings.jsonl")
        self.max_examples = max_examples
        self.counts: Counter[str] = Counter()
        self.example_counts: Counter[str] = Counter()
        self.started = time.monotonic()
        output.parent.mkdir(parents=True, exist_ok=True)
        self.findings = self.findings_path.open("w", encoding="utf-8")

    def compare(
        self,
        feature: str,
        record: dict[str, Any],
        rdkit_value: Any,
        cosmolkit_value: Any,
    ) -> None:
        if rdkit_value == cosmolkit_value:
            self.counts[f"match.{feature}"] += 1
            return
        key = f"mismatch.{feature}"
        self.counts[key] += 1
        if self.example_counts[key] >= self.max_examples:
            return
        self.example_counts[key] += 1
        self.findings.write(
            json.dumps(
                {
                    "feature": feature,
                    "row": record["row"],
                    "chembl_id": record["chembl_id"],
                    "smiles": record["smiles"],
                    "rdkit": rdkit_value,
                    "cosmolkit": cosmolkit_value,
                },
                ensure_ascii=True,
                separators=(",", ":"),
            )
            + "\n"
        )
        self.findings.flush()

    def finish(self, processed: int, profiles: dict[str, int]) -> None:
        self.findings.close()
        temporary = self.output.with_suffix(self.output.suffix + f".tmp-{os.getpid()}")
        temporary.write_text(
            json.dumps(
                {
                    "processed": processed,
                    "elapsed_seconds": time.monotonic() - self.started,
                    "profiles": profiles,
                    "counts": dict(sorted(self.counts.items())),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        temporary.replace(self.output)


def outcome(branch: dict[str, Any]) -> dict[str, Any]:
    error = branch.get("error")
    return {
        "ok": branch.get("ok"),
        "error_type": error.get("type") if isinstance(error, dict) else None,
    }


def compare_record(
    audit: Audit,
    record: dict[str, Any],
    rdkit: dict[str, Any],
    cosmolkit: dict[str, Any],
    branch_names: list[str],
    branch_fields: list[str],
) -> None:
    rdkit_parse = bool(rdkit["parse"]["ok"])
    cosmolkit_parse = bool(cosmolkit["parse"]["ok"])
    audit.compare("tautomer.parse", record, rdkit_parse, cosmolkit_parse)
    if not rdkit_parse or not cosmolkit_parse:
        return
    for name in branch_names:
        rdkit_branch = rdkit["branches"][name]
        cosmolkit_branch = cosmolkit["branches"][name]["result"]
        audit.compare(
            f"tautomer.{name}.enumeration_outcome",
            record,
            outcome(rdkit_branch),
            outcome(cosmolkit_branch),
        )
        if not rdkit_branch["ok"] or not cosmolkit_branch["ok"]:
            continue
        for field in branch_fields:
            audit.compare(
                f"tautomer.{name}.{field}",
                record,
                rdkit_branch[field],
                cosmolkit_branch[field],
            )


def compare_branch(
    audit: Audit,
    record: dict[str, Any],
    name: str,
    rdkit_branch: dict[str, Any],
    cosmolkit_branch: dict[str, Any],
    branch_fields: list[str],
) -> None:
    audit.compare(
        f"tautomer.{name}.enumeration_outcome",
        record,
        outcome(rdkit_branch),
        outcome(cosmolkit_branch),
    )
    if not rdkit_branch["ok"] or not cosmolkit_branch["ok"]:
        return
    for field in branch_fields:
        audit.compare(
            f"tautomer.{name}.{field}",
            record,
            rdkit_branch[field],
            cosmolkit_branch[field],
        )


def read_wire_message(
    stream: Any,
    *,
    kind: str,
    row: int,
    branch: str | None = None,
) -> dict[str, Any]:
    line = stream.readline()
    if not line:
        raise ValueError(f"Rust oracle ended before row {row} {kind}")
    message = json.loads(line)
    if message.get("kind") != kind or message.get("row") != row:
        raise ValueError(f"Rust oracle wire order mismatch at row {row} {kind}")
    if branch is not None and message.get("name") != branch:
        raise ValueError(
            f"Rust oracle branch order mismatch at row {row}: "
            f"expected {branch}, got {message.get('name')}"
        )
    return message


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mode", choices=["tautomer"], required=True)
    parser.add_argument("--tautomer-profile", type=Path, required=True)
    parser.add_argument("--cosmolkit-oracle", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--max-examples", type=int, default=12)
    args = parser.parse_args()

    oracle.assert_rdkit_version()
    profile = oracle.load_profile(args.tautomer_profile)
    branch_names = [str(name) for name in profile["chembl_branches"]]
    branch_fields = [
        str(name)
        for name in profile["comparison_fields"]
        if name not in {"parse", "enumeration_outcome"}
    ]
    branches_by_name = {
        str(branch["name"]): branch for branch in profile["branches"]
    }
    command = [
        str(args.cosmolkit_oracle),
        "--input",
        str(args.input),
        "--output",
        "-",
        "--profile",
        str(args.tautomer_profile),
    ]
    if args.limit is not None:
        command.extend(("--limit", str(args.limit)))
    process = subprocess.Popen(  # noqa: S603
        command,
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
    )
    assert process.stdout is not None
    try:
        audit = Audit(args.output, args.max_examples)
        processed = 0
        with args.input.open(encoding="utf-8") as source:
            records = (line for line in source if line.strip())
            if args.limit is not None:
                records = itertools.islice(records, args.limit)
            for expected_line in records:
                record = json.loads(expected_line)
                reference_input = {
                    **record,
                    "case_id": record.get("chembl_id", f"chembl:{record['row']}"),
                    "sanitize": True,
                    "remove_hs": True,
                    "source": "ChEMBL 37",
                }
                try:
                    rdkit_molecule = oracle.parse_molecule(reference_input)
                except Exception:  # noqa: BLE001
                    rdkit_molecule = None
                parse_message = read_wire_message(
                    process.stdout, kind="parse", row=int(record["row"])
                )
                rdkit_parse = rdkit_molecule is not None
                cosmolkit_parse = bool(parse_message["parse"]["ok"])
                audit.compare(
                    "tautomer.parse", record, rdkit_parse, cosmolkit_parse
                )
                if cosmolkit_parse:
                    for name in branch_names:
                        cosmolkit_branch = read_wire_message(
                            process.stdout,
                            kind="branch",
                            row=int(record["row"]),
                            branch=name,
                        )["result"]
                        if rdkit_parse:
                            assert rdkit_molecule is not None
                            rdkit_branch = oracle.enumerate_branch(
                                rdkit_molecule, branches_by_name[name]
                            )
                            compare_branch(
                                audit,
                                record,
                                name,
                                rdkit_branch,
                                cosmolkit_branch,
                                branch_fields,
                            )
                processed += 1
        if process.stdout.readline():
            raise ValueError("Rust oracle produced records beyond the requested input")
        returncode = process.wait()
        if returncode != 0:
            raise subprocess.CalledProcessError(returncode, command)
        profiles = {"tautomer_parse": 1}
        profiles.update(
            {f"tautomer_{name}": 1 + len(branch_fields) for name in branch_names}
        )
        audit.finish(processed, profiles)
    finally:
        process.stdout.close()
        if process.poll() is None:
            process.terminate()
            process.wait()


if __name__ == "__main__":
    main()
