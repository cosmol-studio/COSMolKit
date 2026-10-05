#!/usr/bin/env python3
"""Prepare the fixed long-conjugated regression using the existing RDKit oracle."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from _tautomer_oracle import assert_rdkit_version, build_record, load_profile


def generate(input_path: Path, output: Path) -> None:
    assert_rdkit_version()
    fixture = json.loads(input_path.read_text(encoding="utf-8"))
    profile = load_profile()
    if (
        fixture["schema_version"] != 1
        or fixture["reference"]["version"] != profile["rdkit_version"]
        or fixture["reference"]["source_revision"] != profile["source_revision"]
        or fixture["branches"] != profile["branches"][:2]
        or len(fixture["cases"]) != 1
    ):
        raise ValueError("long-conjugated fixture/pinned profile mismatch")
    # No CK output participates in preparation; preserve both original branches.
    rows = [build_record(case, profile, ["default", "v1"]) for case in fixture["cases"]]
    output.write_text(
        "".join(json.dumps(row, ensure_ascii=True, sort_keys=True, separators=(",", ":")) + "\n" for row in rows),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    generate(args.input, args.output)


if __name__ == "__main__":
    main()
